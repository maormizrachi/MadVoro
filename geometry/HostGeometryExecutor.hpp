#ifndef MADVORO_GEOMETRY_HOST_GEOMETRY_EXECUTOR_HPP
#define MADVORO_GEOMETRY_HOST_GEOMETRY_EXECUTOR_HPP

#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

#ifdef MADVORO_WITH_KOKKOS
#include <Kokkos_Core.hpp>
#include "../portability/KokkosInitializer.hpp"
#endif

#include "FlatMesh.hpp"
#include "FlatVoronoi.hpp"
#include "kernels/CellMoments.hpp"
#include "kernels/Circumcenter.hpp"
#include "kernels/FaceBuild.hpp"
#include "../utils/Predicates3D.hpp"

namespace MadVoro
{
namespace Geometry
{

#ifdef MADVORO_WITH_KOKKOS
namespace executor_detail
{
    using HostPolicy = Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic>,
                                           Kokkos::IndexType<std::int64_t>>;
} // namespace executor_detail
#endif

/*! \brief Runs f(i) for i in [0, n) on the host: on the threads of Kokkos's host
    execution space when MadVoro is built with Kokkos, one after the other otherwise.
    f must not throw. */
template<typename F>
void HostParallelFor(flat_index_t n, const F &f)
{
#ifdef MADVORO_WITH_KOKKOS
    KokkosInitializer::Instance();
    Kokkos::parallel_for("MadVoro::HostParallelFor", executor_detail::HostPolicy(0, static_cast<std::int64_t>(n)),
                         [&f](const std::int64_t i) { f(static_cast<flat_index_t>(i)); });
    Kokkos::DefaultHostExecutionSpace().fence();
#else
    for(flat_index_t i = 0; i < n; ++i)
    {
        f(i);
    }
#endif
}

//! \brief A face-building failure, at the earliest point (in build order) where it occurred.
struct FaceBuildReport
{
    FaceBuildResult result;
    flat_index_t point = FLAT_INVALID_INDEX;
};

namespace executor_detail
{
    // Per-thread scratch, sized so that BuildPointFaces cannot overflow it: a ring
    // around one edge is at most the point's tetrahedra plus one, and a point has
    // at most three neighbors per tetrahedron.
    struct HostFaceScratch
    {
        std::vector<flat_index_t> fan, face, backup, flip, accepted;
        std::vector<FlatPoint> points;
        std::vector<double> diffs, areas;
        FaceScratch view;

        explicit HostFaceScratch(flat_index_t maxPointTetras)
        {
            const flat_index_t capacity = maxPointTetras + 2;
            fan.resize(capacity);
            face.resize(capacity);
            backup.resize(capacity);
            flip.resize(capacity);
            points.resize(capacity);
            diffs.resize(capacity);
            areas.resize(capacity);
            accepted.resize(3 * maxPointTetras + 3);
            view.fan = fan.data();
            view.face = face.data();
            view.backup = backup.data();
            view.flip = flip.data();
            view.points = points.data();
            view.diffs = diffs.data();
            view.areas = areas.data();
            view.capacity = capacity;
            view.accepted = accepted.data();
            view.acceptedCapacity = accepted.size();
        }
    };

    // The faces one thread built, appended point after point.
    struct ThreadFaces
    {
        std::vector<flat_index_t> neighbor1, vertexCount, vertices;
        std::vector<double> area;
        std::vector<FlatPoint> cm;
    };

    // Faces of the points a thread processes, in its own buffer; the neighbor0 of
    // each face is its point, so it is not stored.
    struct ThreadFacesSink
    {
        ThreadFaces *out;

        void Add(flat_index_t, flat_index_t other, const flat_index_t *faceVertices, flat_index_t count,
                 double area, const FlatPoint &cm)
        {
            this->out->neighbor1.push_back(other);
            this->out->vertexCount.push_back(count);
            this->out->vertices.insert(this->out->vertices.end(), faceVertices, faceVertices + count);
            this->out->area.push_back(area);
            this->out->cm.push_back(cm);
        }
    };

    // Where a point's faces are in its thread's buffer.
    struct PointFaces
    {
        int thread = 0;
        flat_index_t firstFace = 0;
        flat_index_t faceCount = 0;
        flat_index_t firstVertex = 0;
        flat_index_t vertexCount = 0;
    };

    inline flat_index_t MaxPointTetras(const FlatMeshView &mesh)
    {
        flat_index_t result = 0;
        for(flat_index_t p = 0; p < mesh.ownedCount; ++p)
        {
            result = (mesh.PointTetraCount(p) > result) ? mesh.PointTetraCount(p) : result;
        }
        return result;
    }
} // namespace executor_detail

/*! \brief Builds every face of the owned points into `out`, on the host. Points are
    processed in `order`, in parallel when MadVoro is built with Kokkos; each thread
    appends the faces of its points to its own buffer, and the buffers are then copied into place in `order`. Faces are
    therefore numbered as if built one point after the other in `order`, and each cell
    lists its faces in increasing face order.
    `out` must already hold the circumcenters (out.vertices).
    \return The failure at the earliest point in `order`, if any (nothing is filled then). */
inline FaceBuildReport BuildFacesHost(const FlatMeshView &mesh, const flat_index_t *order,
                                      const std::uint8_t *looseAreaCut, FlatVoronoi &out)
{
    const flat_index_t pointCount = mesh.ownedCount;
    const flat_index_t maxPointTetras = executor_detail::MaxPointTetras(mesh);
    const FlatPoint *centers = out.vertices.data();
    std::vector<FaceBuildResult> results(pointCount);
    std::vector<executor_detail::PointFaces> pointFaces(pointCount);
    // One face buffer and scratch per concurrently running thread (a Kokkos unique token)
#ifdef MADVORO_WITH_KOKKOS
    KokkosInitializer::Instance();
    const Kokkos::Experimental::UniqueToken<Kokkos::DefaultHostExecutionSpace> tokens;
    const int slotCount = tokens.size();
#else
    const int slotCount = 1;
#endif
    std::vector<executor_detail::ThreadFaces> buffers(slotCount);
    std::vector<executor_detail::HostFaceScratch> scratches;
    scratches.reserve(slotCount);
    for(int slot = 0; slot < slotCount; ++slot)
        scratches.emplace_back(maxPointTetras);
    const auto buildPoint = [&](flat_index_t i, int slot) {
        executor_detail::ThreadFaces &buffer = buffers[slot];
        executor_detail::ThreadFacesSink sink{&buffer};
        executor_detail::PointFaces &record = pointFaces[i];
        record.thread = slot;
        record.firstFace = buffer.area.size();
        record.firstVertex = buffer.vertices.size();
        results[i] = BuildPointFaces(order[i], mesh, centers, looseAreaCut, scratches[slot].view, sink);
        record.faceCount = buffer.area.size() - record.firstFace;
        record.vertexCount = buffer.vertices.size() - record.firstVertex;
    };
#ifdef MADVORO_WITH_KOKKOS
    Kokkos::parallel_for("MadVoro::BuildFacesHost", executor_detail::HostPolicy(0, static_cast<std::int64_t>(pointCount)),
                         [&](const std::int64_t i) {
                             const int slot = tokens.acquire();
                             buildPoint(static_cast<flat_index_t>(i), slot);
                             tokens.release(slot);
                         });
    Kokkos::DefaultHostExecutionSpace().fence();
#else
    for(flat_index_t i = 0; i < pointCount; ++i)
        buildPoint(i, 0);
#endif
    FaceBuildReport report;
    for(flat_index_t i = 0; i < pointCount; ++i)
    {
        if(results[i].error != FaceBuildError::None)
        {
            report.result = results[i];
            report.point = order[i];
            return report;
        }
    }

    // Each point's first face and vertex, in `order`
    std::vector<flat_index_t> faceOffsets(pointCount + 1, 0), vertexOffsets(pointCount + 1, 0);
    for(flat_index_t i = 0; i < pointCount; ++i)
    {
        faceOffsets[i + 1] = faceOffsets[i] + pointFaces[i].faceCount;
        vertexOffsets[i + 1] = vertexOffsets[i] + pointFaces[i].vertexCount;
    }
    const flat_index_t faceCount = faceOffsets[pointCount];
    out.ResizeFaces(faceCount);
    out.faceVertices.assign(vertexOffsets[pointCount], 0);
    out.faceVertexOffsets[faceCount] = vertexOffsets[pointCount];
    out.cellOrder.assign(order, order + pointCount);
    FlatVoronoiView view = out.View();
    HostParallelFor(pointCount, [&](flat_index_t i) {
        const executor_detail::PointFaces &record = pointFaces[i];
        const executor_detail::ThreadFaces &buffer = buffers[record.thread];
        flat_index_t vertex = vertexOffsets[i];
        flat_index_t source = record.firstVertex;
        for(flat_index_t k = 0; k < record.faceCount; ++k)
        {
            const flat_index_t face = faceOffsets[i] + k;
            const flat_index_t local = record.firstFace + k;
            const flat_index_t count = buffer.vertexCount[local];
            view.faceNeighbor0[face] = order[i];
            view.faceNeighbor1[face] = buffer.neighbor1[local];
            view.faceArea[face] = buffer.area[local];
            view.faceCM[face] = buffer.cm[local];
            view.faceVertexOffsets[face] = vertex;
            for(flat_index_t v = 0; v < count; ++v)
                view.faceVertices[vertex + v] = buffer.vertices[source + v];
            vertex += count;
            source += count;
        }
    });

    HostParallelFor(faceCount, [&](flat_index_t f) {
        view.faceCM[f] = ProjectFaceCM(view.faceCM[f], mesh.points[view.faceNeighbor0[f]],
                                       mesh.points[view.faceNeighbor1[f]]);
    });

    // Each owned cell's faces, in increasing face order
    std::vector<flat_index_t> &cellOffsets = out.cellFaceOffsets;
    cellOffsets.assign(pointCount + 1, 0);
    for(flat_index_t f = 0; f < faceCount; ++f)
    {
        ++cellOffsets[view.faceNeighbor0[f] + 1];
        if(view.faceNeighbor1[f] < pointCount)
            ++cellOffsets[view.faceNeighbor1[f] + 1];
    }
    for(flat_index_t c = 0; c < pointCount; ++c)
        cellOffsets[c + 1] += cellOffsets[c];
    out.cellFaces.assign(cellOffsets[pointCount], 0);
    std::vector<flat_index_t> next(cellOffsets.begin(), cellOffsets.end() - 1);
    for(flat_index_t f = 0; f < faceCount; ++f)
    {
        out.cellFaces[next[view.faceNeighbor0[f]]++] = f;
        if(view.faceNeighbor1[f] < pointCount)
            out.cellFaces[next[view.faceNeighbor1[f]]++] = f;
    }
    return report;
}

/*! \brief Recomputes, with the exact orient3d, the cells whose status is NeedsHost
    (their volumes could not be certified by the static filter). */
inline void ResolveCellMomentsOnHost(const FlatMeshView &mesh, FlatVoronoi &out)
{
    const ConstFlatVoronoiView view = static_cast<const FlatVoronoi &>(out).View();
    const auto exactVolume = [](const FlatPoint &p0, const FlatPoint &p1, const FlatPoint &p2, const FlatPoint &p3,
                                double &volume) {
        volume = std::abs(orient3d(std::array<FlatPoint, 4>{{p0, p1, p2, p3}})) / 6.0;
        return GeometryStatus::Ok;
    };
    for(flat_index_t c = 0; c < view.cellCount; ++c)
    {
        if(out.cellStatus[c] == GeometryStatus::NeedsHost)
        {
            out.cellStatus[c] = CellMoments(c, mesh.points[c], view, exactVolume, out.cellVolume[c], out.cellCM[c]);
        }
    }
}

/*! \brief Volume and centroid of every owned cell into out.cellVolume and out.cellCM,
    on the host. Cells whose volumes the static filter could not certify are
    recomputed with the exact orient3d. */
inline void CellMomentsHost(const FlatMeshView &mesh, FlatVoronoi &out)
{
    const ConstFlatVoronoiView view = static_cast<const FlatVoronoi &>(out).View();
    // In build order, so that neighboring cells (and the circumcenters they share) are visited together
    const bool ordered = (out.cellOrder.size() == view.cellCount);
    HostParallelFor(view.cellCount, [&](flat_index_t i) {
        const flat_index_t c = ordered ? out.cellOrder[i] : i;
        out.cellStatus[c] = CellMoments(c, mesh.points[c], view, CertifiedTetraVolume{}, out.cellVolume[c], out.cellCM[c]);
    });
    ResolveCellMomentsOnHost(mesh, out);
}

/*! \brief Runs the geometry kernels on the host: on Kokkos's host threads when MadVoro
    is built with Kokkos, one after the other otherwise.
    The calls of one build are, in order: SetMesh, Circumcenters, SetVertices,
    BuildFaces, CellMoments. Inputs and results live in the host FlatMesh and
    FlatVoronoi. */
class HostGeometryExecutor
{
public:
    //! \brief Makes `mesh` the input of the following calls; it must outlive them.
    void SetMesh(const FlatMesh &mesh)
    {
        this->mesh_ = &mesh;
    }

    /*! \brief Circumcenters of the tetrahedra with compute[t] != 0, into out.vertices,
        out.tetraRadius and out.tetraStatus. The other tetrahedra keep their vertices. */
    void Circumcenters(const std::vector<std::uint8_t> &compute, FlatVoronoi &out)
    {
        const FlatMeshView mesh = this->mesh_->View();
        const FlatVoronoiView view = out.View();
        HostParallelFor(compute.size(), [&](flat_index_t t) {
            if(compute[t] != 0)
                CircumcenterOfTetra(t, mesh, view);
        });
    }

    //! \brief Makes out.vertices the face vertices for BuildFaces (nothing to do on the host).
    void SetVertices(const FlatVoronoi &)
    {
    }

    //! \brief See BuildFacesHost.
    FaceBuildReport BuildFaces(const std::vector<flat_index_t> &order, const std::vector<std::uint8_t> &looseAreaCut,
                               FlatVoronoi &out)
    {
        return BuildFacesHost(this->mesh_->View(), order.data(), looseAreaCut.data(), out);
    }

    //! \brief Volumes and centroids of the cells of the last BuildFaces, into `out`.
    void CellMoments(FlatVoronoi &out)
    {
        CellMomentsHost(this->mesh_->View(), out);
    }

private:
    const FlatMesh *mesh_ = nullptr;
};

} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_EXECUTOR_HPP
