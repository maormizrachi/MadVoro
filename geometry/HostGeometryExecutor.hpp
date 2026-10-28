#ifndef MADVORO_GEOMETRY_HOST_GEOMETRY_EXECUTOR_HPP
#define MADVORO_GEOMETRY_HOST_GEOMETRY_EXECUTOR_HPP

#include <array>
#include <cmath>
#include <cstdint>
#include <optional>
#include <vector>

#ifdef MADVORO_WITH_GPU
#include <Kokkos_Core.hpp>
#include "../portability/KokkosInitializer.hpp"
#endif

#include "FlatMesh.hpp"
#include "FlatVoronoi.hpp"
#include "RawFaces.hpp"
#include "kernels/CellMoments.hpp"
#include "kernels/Circumcenter.hpp"
#include "kernels/FaceBuild.hpp"
#include "../utils/Predicates3D.hpp"

namespace MadVoro
{
namespace Geometry
{

#ifdef MADVORO_WITH_GPU
namespace executor_detail
{
    using HostPolicy = Kokkos::RangePolicy<Kokkos::DefaultHostExecutionSpace, Kokkos::Schedule<Kokkos::Dynamic>, Kokkos::IndexType<std::int64_t>>;
} // namespace executor_detail
#endif

namespace executor_detail
{
    inline bool &SerialScopeFlag()
    {
        thread_local bool serial = false;
        return serial;
    }
} // namespace executor_detail

//! \brief While alive, geometry work on this thread runs serially (Kokkos cannot launch a parallel loop from inside one).
class SerialGeometryScope
{
public:
    SerialGeometryScope() : previous_(executor_detail::SerialScopeFlag())
    {
        executor_detail::SerialScopeFlag() = true;
    }

    ~SerialGeometryScope()
    {
        executor_detail::SerialScopeFlag() = this->previous_;
    }

    SerialGeometryScope(const SerialGeometryScope &) = delete;
    SerialGeometryScope &operator=(const SerialGeometryScope &) = delete;

private:
    bool previous_;
};

inline bool InSerialGeometryScope()
{
    return executor_detail::SerialScopeFlag();
}

//! \brief f must not throw.
template<typename F>
void HostParallelFor(flat_index_t n, const F &f)
{
#ifdef MADVORO_WITH_GPU
    if(not InSerialGeometryScope())
    {
        KokkosInitializer::Instance();
        Kokkos::parallel_for("MadVoro::HostParallelFor", executor_detail::HostPolicy(0, static_cast<std::int64_t>(n)), [&f](const std::int64_t i) { f(static_cast<flat_index_t>(i)); });
        Kokkos::DefaultHostExecutionSpace().fence();
        return;
    }
#endif
    for(flat_index_t i = 0; i < n; ++i)
    {
        f(i);
    }
}

//! \brief The failure at the earliest point in build order.
struct FaceBuildReport
{
    FaceBuildResult result;
    flat_index_t point = FLAT_INVALID_INDEX;
};

namespace executor_detail
{
    // Sized so that BuildPointFaces cannot overflow it.
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

    struct ThreadFaces
    {
        std::vector<flat_index_t> other, vertexCount, vertices;
    };

    // a face's cell is its point, so it is not stored.
    struct ThreadFacesSink
    {
        ThreadFaces *out;

        void Add(flat_index_t, flat_index_t other, const flat_index_t *loop, flat_index_t count)
        {
            this->out->other.push_back(other);
            this->out->vertexCount.push_back(count);
            this->out->vertices.insert(this->out->vertices.end(), loop, loop + count);
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

//! \brief The raw faces of the points in `order` (each point's faces together, the points in that order)
inline FaceBuildReport BuildRawFacesHost(const FlatMeshView &mesh, const flat_index_t *order, HostRawFaces &raw)
{
    const flat_index_t pointCount = mesh.ownedCount;
    const flat_index_t maxPointTetras = executor_detail::MaxPointTetras(mesh);
    std::vector<FaceBuildResult> results(pointCount);
    std::vector<executor_detail::PointFaces> pointFaces(pointCount);
    // One face buffer and scratch per concurrently running thread (a Kokkos unique token)
#ifdef MADVORO_WITH_GPU
    const bool parallel = not InSerialGeometryScope();
    std::optional<Kokkos::Experimental::UniqueToken<Kokkos::DefaultHostExecutionSpace>> tokens;
    if(parallel)
    {
        KokkosInitializer::Instance();
        tokens.emplace();
    }
    const int slotCount = parallel ? tokens->size() : 1;
#else
    const bool parallel = false;
    const int slotCount = 1;
#endif
    std::vector<executor_detail::ThreadFaces> buffers(slotCount);
    std::vector<executor_detail::HostFaceScratch> scratches;
    scratches.reserve(slotCount);
    for(int slot = 0; slot < slotCount; ++slot)
    {
        scratches.emplace_back(maxPointTetras);
    }
    const auto buildPoint = [&](flat_index_t i, int slot)
    {
        executor_detail::ThreadFaces &buffer = buffers[slot];
        executor_detail::ThreadFacesSink sink{&buffer};
        executor_detail::PointFaces &record = pointFaces[i];
        record.thread = slot;
        record.firstFace = buffer.other.size();
        record.firstVertex = buffer.vertices.size();
        results[i] = BuildPointRawFaces(order[i], mesh, scratches[slot].view, sink);
        record.faceCount = buffer.other.size() - record.firstFace;
        record.vertexCount = buffer.vertices.size() - record.firstVertex;
    };
#ifdef MADVORO_WITH_GPU
    if(parallel)
    {
        Kokkos::parallel_for("MadVoro::BuildRawFacesHost", executor_detail::HostPolicy(0, static_cast<std::int64_t>(pointCount)), [&](const std::int64_t i)
        {
            const int slot = tokens->acquire();
            buildPoint(static_cast<flat_index_t>(i), slot);
            tokens->release(slot);
        });
        Kokkos::DefaultHostExecutionSpace().fence();
    }
#endif
    if(not parallel)
    {
        for(flat_index_t i = 0; i < pointCount; ++i)
        {
            buildPoint(i, 0);
        }
    }
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
    raw.faceCell.resize(faceCount);
    raw.faceOther.resize(faceCount);
    raw.faceVertexOffsets.resize(faceCount + 1);
    raw.faceVertices.resize(vertexOffsets[pointCount]);
    raw.faceVertexOffsets[faceCount] = vertexOffsets[pointCount];
    HostParallelFor(pointCount, [&](flat_index_t i)
    {
        const executor_detail::PointFaces &record = pointFaces[i];
        const executor_detail::ThreadFaces &buffer = buffers[record.thread];
        flat_index_t vertex = vertexOffsets[i];
        flat_index_t source = record.firstVertex;
        for(flat_index_t k = 0; k < record.faceCount; ++k)
        {
            const flat_index_t face = faceOffsets[i] + k;
            const flat_index_t count = buffer.vertexCount[record.firstFace + k];
            raw.faceCell[face] = order[i];
            raw.faceOther[face] = buffer.other[record.firstFace + k];
            raw.faceVertexOffsets[face] = vertex;
            for(flat_index_t v = 0; v < count; ++v)
            {
                raw.faceVertices[vertex + v] = buffer.vertices[source + v];
            }
            vertex += count;
            source += count;
        }
    });
    return report;
}

//! \brief The view of raw faces on the host
inline RawFacesView HostRawFacesView(const FlatMesh &mesh, const FlatPoint *vertices, flat_index_t vertexCount, const flat_index_t *order, const HostRawFaces &raw)
{
    RawFacesView view;
    view.points = mesh.points.data();
    view.hostPoints = mesh.points.data();
    view.vertices = vertices;
    view.faceCell = raw.faceCell.data();
    view.faceOther = raw.faceOther.data();
    view.faceVertexOffsets = raw.faceVertexOffsets.data();
    view.faceVertices = raw.faceVertices.data();
    view.order = order;
    view.cellCount = mesh.ownedCount;
    view.pointCount = mesh.points.size();
    view.vertexCount = vertexCount;
    view.faceCount = raw.faceCell.size();
    view.faceVertexCount = raw.faceVertices.size();
    view.onDevice = false;
    return view;
}

//! \brief Recomputes the NeedsHost cells with the exact orient3d.
inline void ResolveCellMomentsOnHost(const FlatMeshView &mesh, FlatVoronoi &out)
{
    const ConstFlatVoronoiView view = static_cast<const FlatVoronoi &>(out).View();
    const auto exactVolume = [](const FlatPoint &p0, const FlatPoint &p1, const FlatPoint &p2, const FlatPoint &p3, double &volume)
    {
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

inline void CellMomentsHost(const FlatMeshView &mesh, FlatVoronoi &out)
{
    const ConstFlatVoronoiView view = static_cast<const FlatVoronoi &>(out).View();
    // In build order, so that neighboring cells (and the circumcenters they share) are visited together
    const bool ordered = (out.cellOrder.size() == view.cellCount);
    HostParallelFor(view.cellCount, [&](flat_index_t i)
    {
        const flat_index_t c = ordered ? out.cellOrder[i] : i;
        out.cellStatus[c] = CellMoments(c, mesh.points[c], view, CertifiedTetraVolume{}, out.cellVolume[c], out.cellCM[c]);
    });
    ResolveCellMomentsOnHost(mesh, out);
}

//! \brief Calls of one build, in order: SetMesh, Circumcenters, SetVertices, BuildRawFaces (whose view stays valid until the next SetMesh).
class HostGeometryExecutor
{
public:
    //! \brief `mesh` must outlive the following calls.
    void SetMesh(const FlatMesh &mesh)
    {
        this->mesh_ = &mesh;
    }

    //! \brief Only the tetrahedra with compute[t] != 0; the others keep their vertices.
    void Circumcenters(const std::vector<std::uint8_t> &compute, FlatVoronoi &out)
    {
        const FlatMeshView mesh = this->mesh_->View();
        const FlatVoronoiView view = out.View();
        HostParallelFor(compute.size(), [&](flat_index_t t)
        {
            if(compute[t] != 0)
            {
                CircumcenterOfTetra(t, mesh, view);
            }
        });
    }

    //! \brief `out` (its vertices) must outlive the raw faces.
    void SetVertices(const FlatVoronoi &out)
    {
        this->vertices_ = &out.vertices;
    }

    //! \brief The raw faces of the mesh's points in `order` (which must outlive them)
    FaceBuildReport BuildRawFaces(const std::vector<flat_index_t> &order, RawFacesView &raw)
    {
        const FaceBuildReport report = BuildRawFacesHost(this->mesh_->View(), order.data(), this->raw_);
        raw = HostRawFacesView(*this->mesh_, this->vertices_->data(), this->vertices_->size(), order.data(), this->raw_);
        return report;
    }

    void Release(void)
    {
        this->raw_.Release();
    }

private:
    const FlatMesh *mesh_ = nullptr;
    const std::vector<FlatPoint> *vertices_ = nullptr;
    HostRawFaces raw_;
};

} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_EXECUTOR_HPP
