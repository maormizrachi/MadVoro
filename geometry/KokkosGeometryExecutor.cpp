#include "KokkosGeometryExecutor.hpp"

#include <string>

#include "kernels/CellMoments.hpp"
#include "kernels/Circumcenter.hpp"
#include "kernels/FaceBuild.hpp"
#include "../portability/KokkosInitializer.hpp"

namespace MadVoro
{
namespace Geometry
{

namespace
{
    using Policy = Kokkos::RangePolicy<Kokkos::IndexType<std::int64_t>>;

    // Fixed scratch of the face kernels: a ring around one Delaunay edge, and the
    // neighbors of one point. A point that does not fit sends the build to the host.
    constexpr int FAN_CAPACITY = 64;
    constexpr int ACCEPTED_CAPACITY = 192;
    using KernelFaceScratch = FixedFaceScratch<FAN_CAPACITY, ACCEPTED_CAPACITY>;

    template<typename T>
    void EnsureSize(Kokkos::View<T *> &view, std::size_t size, const char *label)
    {
        if(view.extent(0) < size)
            view = Kokkos::View<T *>(Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string(label)), size);
    }
} // namespace

void KokkosGeometryExecutor::SetMesh(const FlatMesh &mesh)
{
    KokkosInitializer::Instance();
    this->hostMesh_ = &mesh;
    this->points_.Upload(mesh.points);
    this->tetras_.Upload(mesh.tetras);
    this->pointTetraOffsets_.Upload(mesh.pointTetraOffsets);
    this->pointTetras_.Upload(mesh.pointTetras);
    this->mesh_ = mesh.View();
    this->mesh_.points = this->points_.Data();
    this->mesh_.tetras = this->tetras_.Data();
    this->mesh_.pointTetraOffsets = this->pointTetraOffsets_.Data();
    this->mesh_.pointTetras = this->pointTetras_.Data();
}

FlatVoronoiView KokkosGeometryExecutor::DeviceView(const FlatVoronoi &out) const
{
    FlatVoronoiView view;
    view.cellFaceOffsets = this->cellFaceOffsets_.Data();
    view.cellFaces = this->cellFaces_.Data();
    view.faceVertexOffsets = this->faceVertexOffsets_.Data();
    view.faceVertices = this->faceVertices_.Data();
    view.faceNeighbor0 = this->faceNeighbor0_.Data();
    view.faceNeighbor1 = this->faceNeighbor1_.Data();
    view.faceArea = this->faceArea_.Data();
    view.faceCM = this->faceCM_.Data();
    view.cellVolume = this->cellVolume_.Data();
    view.cellCM = this->cellCM_.Data();
    view.cellStatus = this->cellStatus_.Data();
    view.vertices = this->vertices_.Data();
    view.tetraRadius = this->tetraRadius_.Data();
    view.tetraStatus = this->tetraStatus_.Data();
    view.cellCount = out.CellCount();
    view.faceCount = out.FaceCount();
    view.tetraCount = out.TetraCount();
    return view;
}

void KokkosGeometryExecutor::Circumcenters(const std::vector<std::uint8_t> &compute, FlatVoronoi &out)
{
    this->compute_.Upload(compute);
    // Reused circumcenters are already in out.vertices
    this->vertices_.Upload(out.vertices);
    this->tetraRadius_.Upload(out.tetraRadius);
    this->tetraStatus_.Upload(out.tetraStatus);
    const FlatMeshView mesh = this->mesh_;
    const FlatVoronoiView view = this->DeviceView(out);
    const std::uint8_t *mask = this->compute_.Data();
    Kokkos::parallel_for("MadVoro::Circumcenters", Policy(0, static_cast<std::int64_t>(compute.size())),
                         KOKKOS_LAMBDA(const std::int64_t t) {
                             if(mask[t] != 0)
                                 CircumcenterOfTetra(static_cast<flat_index_t>(t), mesh, view);
                         });
    Kokkos::fence();
    this->vertices_.Download(out.vertices);
    this->tetraRadius_.Download(out.tetraRadius);
    this->tetraStatus_.Download(out.tetraStatus);
}

void KokkosGeometryExecutor::SetVertices(const FlatVoronoi &out)
{
    this->vertices_.Upload(out.vertices);
}

void KokkosGeometryExecutor::LaunchFaceCount(std::size_t pointCount)
{
    const FlatMeshView mesh = this->mesh_;
    const FlatPoint *centers = this->vertices_.Data();
    const std::uint8_t *looseAreaCut = this->looseAreaCut_.Data();
    const flat_index_t *order = this->order_.Data();
    const Kokkos::View<FaceBuildResult *> results = this->results_;
    const Kokkos::View<flat_index_t *> faceOffsets = this->faceOffsets_;
    const Kokkos::View<flat_index_t *> vertexOffsets = this->vertexOffsets_;
    const std::int64_t count = static_cast<std::int64_t>(pointCount);

    // Pass 1: each point's faces and face vertices, at index i + 1 ...
    Kokkos::parallel_for("MadVoro::FaceCount", Policy(0, count), KOKKOS_LAMBDA(const std::int64_t i) {
        KernelFaceScratch scratch;
        FaceCountSink sink;
        results(i) = BuildPointFaces(order[i], mesh, centers, looseAreaCut, scratch.View(), sink);
        faceOffsets(i + 1) = sink.faces;
        vertexOffsets(i + 1) = sink.vertices;
        if(i == 0)
        {
            faceOffsets(0) = 0;
            vertexOffsets(0) = 0;
        }
    });
    // ... turned into each point's first face and vertex
    Kokkos::parallel_scan("MadVoro::FaceOffsets", Policy(0, count),
                          KOKKOS_LAMBDA(const std::int64_t i, flat_index_t &update, const bool final) {
                              update += faceOffsets(i + 1);
                              if(final)
                                  faceOffsets(i + 1) = update;
                          });
    Kokkos::parallel_scan("MadVoro::VertexOffsets", Policy(0, count),
                          KOKKOS_LAMBDA(const std::int64_t i, flat_index_t &update, const bool final) {
                              update += vertexOffsets(i + 1);
                              if(final)
                                  vertexOffsets(i + 1) = update;
                          });
    Kokkos::fence();
}

void KokkosGeometryExecutor::LaunchFaceFill(std::size_t pointCount, const FlatVoronoiView &view)
{
    const FlatMeshView mesh = this->mesh_;
    const FlatPoint *centers = this->vertices_.Data();
    const std::uint8_t *looseAreaCut = this->looseAreaCut_.Data();
    const flat_index_t *order = this->order_.Data();
    const Kokkos::View<flat_index_t *> faceOffsets = this->faceOffsets_;
    const Kokkos::View<flat_index_t *> vertexOffsets = this->vertexOffsets_;
    const FlatVoronoiView out = view;

    // Pass 2: the same faces, written from each point's first face and vertex
    Kokkos::parallel_for("MadVoro::FaceFill", Policy(0, static_cast<std::int64_t>(pointCount)),
                         KOKKOS_LAMBDA(const std::int64_t i) {
                             KernelFaceScratch scratch;
                             FaceFillSink sink{out, faceOffsets(i), vertexOffsets(i)};
                             BuildPointFaces(order[i], mesh, centers, looseAreaCut, scratch.View(), sink);
                         });
    Kokkos::parallel_for("MadVoro::FaceCM", Policy(0, static_cast<std::int64_t>(out.faceCount)),
                         KOKKOS_LAMBDA(const std::int64_t f) {
                             out.faceCM[f] = ProjectFaceCM(out.faceCM[f], mesh.points[out.faceNeighbor0[f]],
                                                           mesh.points[out.faceNeighbor1[f]]);
                         });
    Kokkos::fence();
}

void KokkosGeometryExecutor::LaunchCellFaceLists(std::size_t pointCount, std::size_t faceCount, FlatVoronoi &hostOut)
{
    const flat_index_t cells = pointCount;
    const flat_index_t *neighbor0 = this->faceNeighbor0_.Data();
    const flat_index_t *neighbor1 = this->faceNeighbor1_.Data();
    const std::int64_t faces = static_cast<std::int64_t>(faceCount);

    // Faces per owned cell, then offsets
    Kokkos::View<flat_index_t *> counts(std::string("MadVoro::CellFaceCounts"), cells);
    Kokkos::parallel_for("MadVoro::CellFaceCount", Policy(0, faces), KOKKOS_LAMBDA(const std::int64_t f) {
        Kokkos::atomic_inc(&counts(neighbor0[f]));
        if(neighbor1[f] < cells)
            Kokkos::atomic_inc(&counts(neighbor1[f]));
    });
    hostOut.cellFaceOffsets.assign(cells + 1, 0);
    this->cellFaceOffsets_.Upload(hostOut.cellFaceOffsets);
    flat_index_t *offsets = this->cellFaceOffsets_.Data();
    flat_index_t total = 0;
    Kokkos::parallel_scan("MadVoro::CellFaceOffsets", Policy(0, static_cast<std::int64_t>(cells)),
                          KOKKOS_LAMBDA(const std::int64_t i, flat_index_t &update, const bool final) {
                              update += counts(i);
                              if(final)
                                  offsets[i + 1] = update;
                          }, total);

    // Fill, then sort each cell's faces into increasing face order (the fill order is arbitrary)
    hostOut.cellFaces.assign(total, 0);
    this->cellFaces_.Output(hostOut.cellFaces);
    flat_index_t *cellFaces = this->cellFaces_.Data();
    Kokkos::View<flat_index_t *> cursor(Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string("MadVoro::CellFaceCursor")), cells);
    Kokkos::parallel_for("MadVoro::CellFaceCursor", Policy(0, static_cast<std::int64_t>(cells)),
                         KOKKOS_LAMBDA(const std::int64_t c) { cursor(c) = offsets[c]; });
    Kokkos::parallel_for("MadVoro::CellFaceFill", Policy(0, faces), KOKKOS_LAMBDA(const std::int64_t f) {
        cellFaces[Kokkos::atomic_fetch_add(&cursor(neighbor0[f]), flat_index_t(1))] = static_cast<flat_index_t>(f);
        if(neighbor1[f] < cells)
            cellFaces[Kokkos::atomic_fetch_add(&cursor(neighbor1[f]), flat_index_t(1))] = static_cast<flat_index_t>(f);
    });
    Kokkos::parallel_for("MadVoro::CellFaceSort", Policy(0, static_cast<std::int64_t>(cells)),
                         KOKKOS_LAMBDA(const std::int64_t c) {
                             const flat_index_t begin = offsets[c];
                             const flat_index_t end = offsets[c + 1];
                             for(flat_index_t k = begin + 1; k < end; ++k)
                             {
                                 const flat_index_t face = cellFaces[k];
                                 flat_index_t j = k;
                                 while(j > begin && cellFaces[j - 1] > face)
                                 {
                                     cellFaces[j] = cellFaces[j - 1];
                                     --j;
                                 }
                                 cellFaces[j] = face;
                             }
                         });
    Kokkos::fence();
    this->cellFaceOffsets_.Download(hostOut.cellFaceOffsets);
    this->cellFaces_.Download(hostOut.cellFaces);
}

FaceBuildReport KokkosGeometryExecutor::BuildFaces(const std::vector<flat_index_t> &order,
                                                   const std::vector<std::uint8_t> &looseAreaCut, FlatVoronoi &out)
{
    const std::size_t pointCount = this->mesh_.ownedCount;
    this->order_.Upload(order);
    this->looseAreaCut_.Upload(looseAreaCut);
    EnsureSize(this->results_, pointCount, "MadVoro::FaceBuildResults");
    EnsureSize(this->faceOffsets_, pointCount + 1, "MadVoro::FaceOffsets");
    EnsureSize(this->vertexOffsets_, pointCount + 1, "MadVoro::VertexOffsets");
    if(pointCount == 0)
    {
        out.ResizeFaces(0);
        out.faceVertices.clear();
        out.cellFaces.clear();
        out.cellFaceOffsets.assign(1, 0);
        out.cellOrder.clear();
        return FaceBuildReport();
    }

    this->LaunchFaceCount(pointCount);

    // The earliest failure in build order decides, as on the host. A point that did
    // not fit the fixed scratch is built by the host instead (with all others).
    const auto results = Kokkos::create_mirror_view_and_copy(
        Kokkos::HostSpace(), Kokkos::subview(this->results_, std::make_pair(std::size_t(0), pointCount)));
    for(std::size_t i = 0; i < pointCount; ++i)
    {
        if(results(i).error == FaceBuildError::None)
            continue;
        if(results(i).error == FaceBuildError::ScratchOverflow)
        {
            const FaceBuildReport report = BuildFacesHost(this->hostMesh_->View(), order.data(), looseAreaCut.data(), out);
            if(report.result.error == FaceBuildError::None)
            {
                this->faceNeighbor0_.Upload(out.faceNeighbor0);
                this->faceNeighbor1_.Upload(out.faceNeighbor1);
                this->faceArea_.Upload(out.faceArea);
                this->faceCM_.Upload(out.faceCM);
                this->faceVertexOffsets_.Upload(out.faceVertexOffsets);
                this->faceVertices_.Upload(out.faceVertices);
                this->cellFaceOffsets_.Upload(out.cellFaceOffsets);
                this->cellFaces_.Upload(out.cellFaces);
            }
            return report;
        }
        FaceBuildReport report;
        report.result = results(i);
        report.point = order[i];
        return report;
    }

    flat_index_t faceCount = 0, vertexCount = 0;
    Kokkos::deep_copy(faceCount, Kokkos::subview(this->faceOffsets_, pointCount));
    Kokkos::deep_copy(vertexCount, Kokkos::subview(this->vertexOffsets_, pointCount));
    out.ResizeFaces(faceCount);
    out.faceVertices.assign(vertexCount, 0);
    out.faceVertexOffsets[faceCount] = vertexCount;
    out.cellOrder.assign(order.begin(), order.end());
    this->faceNeighbor0_.Output(out.faceNeighbor0);
    this->faceNeighbor1_.Output(out.faceNeighbor1);
    this->faceArea_.Output(out.faceArea);
    this->faceCM_.Output(out.faceCM);
    this->faceVertices_.Output(out.faceVertices);
    this->faceVertexOffsets_.Upload(out.faceVertexOffsets);

    this->LaunchFaceFill(pointCount, this->DeviceView(out));
    this->LaunchCellFaceLists(pointCount, faceCount, out);

    this->faceNeighbor0_.Download(out.faceNeighbor0);
    this->faceNeighbor1_.Download(out.faceNeighbor1);
    this->faceArea_.Download(out.faceArea);
    this->faceCM_.Download(out.faceCM);
    this->faceVertices_.Download(out.faceVertices);
    this->faceVertexOffsets_.Download(out.faceVertexOffsets);
    return FaceBuildReport();
}

void KokkosGeometryExecutor::CellMoments(FlatVoronoi &out)
{
    this->cellVolume_.Output(out.cellVolume);
    this->cellCM_.Output(out.cellCM);
    this->cellStatus_.Output(out.cellStatus);
    this->cellOrder_.Upload(out.cellOrder);
    const ConstFlatVoronoiView view = this->DeviceView(out);
    const FlatMeshView mesh = this->mesh_;
    const flat_index_t *cellOrder = this->cellOrder_.Data();
    const bool ordered = (out.cellOrder.size() == out.CellCount());
    double *volume = this->cellVolume_.Data();
    FlatPoint *cm = this->cellCM_.Data();
    GeometryStatus *status = this->cellStatus_.Data();
    Kokkos::parallel_for("MadVoro::CellMoments", Policy(0, static_cast<std::int64_t>(out.CellCount())),
                         KOKKOS_LAMBDA(const std::int64_t i) {
                             const flat_index_t c = ordered ? cellOrder[i] : static_cast<flat_index_t>(i);
                             status[c] = Geometry::CellMoments(c, mesh.points[c], view, CertifiedTetraVolume{},
                                                               volume[c], cm[c]);
                         });
    Kokkos::fence();
    this->cellVolume_.Download(out.cellVolume);
    this->cellCM_.Download(out.cellCM);
    this->cellStatus_.Download(out.cellStatus);
    ResolveCellMomentsOnHost(this->hostMesh_->View(), out);
}

} // namespace Geometry
} // namespace MadVoro
