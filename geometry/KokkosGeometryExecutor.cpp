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

    // A point that does not fit the fixed scratch sends the build to the host.
    constexpr int FAN_CAPACITY = 64;
    constexpr int ACCEPTED_CAPACITY = 192;
    using KernelFaceScratch = FixedFaceScratch<FAN_CAPACITY, ACCEPTED_CAPACITY>;

    template<typename T>
    void EnsureSize(Kokkos::View<T *> &view, std::size_t size, const char *label)
    {
        if(view.extent(0) < size)
        {
            view = Kokkos::View<T *>(Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string(label)), size);
        }
    }
} // namespace

void KokkosGeometryExecutor::SetMesh(const FlatMesh &mesh)
{
    this->host_.SetMesh(mesh);
    if(InSerialGeometryScope())
    {
        return;
    }
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
    if(InSerialGeometryScope())
    {
        this->host_.Circumcenters(compute, out);
        return;
    }
    this->compute_.Upload(compute);
    // Reused circumcenters are already in out.vertices
    this->vertices_.Upload(out.vertices);
    this->tetraRadius_.Upload(out.tetraRadius);
    this->tetraStatus_.Upload(out.tetraStatus);
    const FlatMeshView mesh = this->mesh_;
    const FlatVoronoiView view = this->DeviceView(out);
    const std::uint8_t *mask = this->compute_.Data();
    Kokkos::parallel_for("MadVoro::Circumcenters", Policy(0, static_cast<std::int64_t>(compute.size())), KOKKOS_LAMBDA(const std::int64_t t)
    {
        if(mask[t] != 0)
        {
            CircumcenterOfTetra(static_cast<flat_index_t>(t), mesh, view);
        }
    });
    Kokkos::fence();
    this->vertices_.Download(out.vertices);
    this->tetraRadius_.Download(out.tetraRadius);
    this->tetraStatus_.Download(out.tetraStatus);
}

void KokkosGeometryExecutor::SetVertices(const FlatVoronoi &out)
{
    this->hostVertices_ = &out.vertices;
    if(InSerialGeometryScope())
    {
        this->host_.SetVertices(out);
        return;
    }
    this->vertices_.Upload(out.vertices);
}

void KokkosGeometryExecutor::LaunchFaceCount(std::size_t pointCount)
{
    const FlatMeshView mesh = this->mesh_;
    const flat_index_t *order = this->order_.Data();
    const Kokkos::View<FaceBuildResult *> results = this->results_;
    const Kokkos::View<flat_index_t *> faceOffsets = this->faceOffsets_;
    const Kokkos::View<flat_index_t *> vertexOffsets = this->vertexOffsets_;
    const std::int64_t count = static_cast<std::int64_t>(pointCount);

    // Pass 1: each point's raw faces and their vertices, at index i + 1 ...
    Kokkos::parallel_for("MadVoro::RawFaceCount", Policy(0, count), KOKKOS_LAMBDA(const std::int64_t i)
    {
        KernelFaceScratch scratch;
        RawFaceCountSink sink;
        results(i) = BuildPointRawFaces(order[i], mesh, scratch.View(), sink);
        faceOffsets(i + 1) = sink.faces;
        vertexOffsets(i + 1) = sink.vertices;
        if(i == 0)
        {
            faceOffsets(0) = 0;
            vertexOffsets(0) = 0;
        }
    });
    // ... turned into each point's first face and vertex
    Kokkos::parallel_scan("MadVoro::RawFaceOffsets", Policy(0, count), KOKKOS_LAMBDA(const std::int64_t i, flat_index_t &update, const bool final)
    {
        update += faceOffsets(i + 1);
        if(final)
        {
            faceOffsets(i + 1) = update;
        }
    });
    Kokkos::parallel_scan("MadVoro::RawVertexOffsets", Policy(0, count), KOKKOS_LAMBDA(const std::int64_t i, flat_index_t &update, const bool final)
    {
        update += vertexOffsets(i + 1);
        if(final)
        {
            vertexOffsets(i + 1) = update;
        }
    });
    Kokkos::fence();
}

void KokkosGeometryExecutor::LaunchFaceFill(std::size_t pointCount)
{
    const FlatMeshView mesh = this->mesh_;
    const flat_index_t *order = this->order_.Data();
    const Kokkos::View<flat_index_t *> faceOffsets = this->faceOffsets_;
    const Kokkos::View<flat_index_t *> vertexOffsets = this->vertexOffsets_;
    flat_index_t *rawCell = this->rawCell_.data(), *rawOther = this->rawOther_.data(), *rawVertexOffsets = this->rawVertexOffsets_.data(), *rawVertices = this->rawVertices_.data();

    // Pass 2: the same raw faces, written from each point's first face and vertex
    Kokkos::parallel_for("MadVoro::RawFaceFill", Policy(0, static_cast<std::int64_t>(pointCount)), KOKKOS_LAMBDA(const std::int64_t i)
    {
        KernelFaceScratch scratch;
        RawFaceFillSink sink{rawCell, rawOther, rawVertexOffsets, rawVertices, faceOffsets(i), vertexOffsets(i)};
        BuildPointRawFaces(order[i], mesh, scratch.View(), sink);
    });
    Kokkos::fence();
}

void device_detail::BuildCellFaceLists(flat_index_t cells, flat_index_t faceCount, const flat_index_t *neighbor0, const flat_index_t *neighbor1, DeviceArray<flat_index_t> &cellFaceOffsets, DeviceArray<flat_index_t> &cellFaces, FlatVoronoi &hostOut)
{
    const std::int64_t faces = static_cast<std::int64_t>(faceCount);

    // Faces per owned cell, then offsets
    Kokkos::View<flat_index_t *> counts(std::string("MadVoro::CellFaceCounts"), cells);
    Kokkos::parallel_for("MadVoro::CellFaceCount", Policy(0, faces), KOKKOS_LAMBDA(const std::int64_t f)
    {
        Kokkos::atomic_inc(&counts(neighbor0[f]));
        if(neighbor1[f] < cells)
        {
            Kokkos::atomic_inc(&counts(neighbor1[f]));
        }
    });
    hostOut.cellFaceOffsets.assign(cells + 1, 0);
    cellFaceOffsets.Upload(hostOut.cellFaceOffsets);
    flat_index_t *offsets = cellFaceOffsets.Data();
    flat_index_t total = 0;
    Kokkos::parallel_scan("MadVoro::CellFaceOffsets", Policy(0, static_cast<std::int64_t>(cells)), KOKKOS_LAMBDA(const std::int64_t i, flat_index_t &update, const bool final)
    {
        update += counts(i);
        if(final)
        {
            offsets[i + 1] = update;
        }
    }, total);

    // Fill, then sort each cell's faces into increasing face order (the fill order is arbitrary)
    hostOut.cellFaces.assign(total, 0);
    cellFaces.Output(hostOut.cellFaces);
    flat_index_t *faceList = cellFaces.Data();
    Kokkos::View<flat_index_t *> cursor(Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string("MadVoro::CellFaceCursor")), cells);
    Kokkos::parallel_for("MadVoro::CellFaceCursor", Policy(0, static_cast<std::int64_t>(cells)), KOKKOS_LAMBDA(const std::int64_t c) { cursor(c) = offsets[c]; });
    Kokkos::parallel_for("MadVoro::CellFaceFill", Policy(0, faces), KOKKOS_LAMBDA(const std::int64_t f)
    {
        faceList[Kokkos::atomic_fetch_add(&cursor(neighbor0[f]), flat_index_t(1))] = static_cast<flat_index_t>(f);
        if(neighbor1[f] < cells)
        {
            faceList[Kokkos::atomic_fetch_add(&cursor(neighbor1[f]), flat_index_t(1))] = static_cast<flat_index_t>(f);
        }
    });
    Kokkos::parallel_for("MadVoro::CellFaceSort", Policy(0, static_cast<std::int64_t>(cells)), KOKKOS_LAMBDA(const std::int64_t c)
    {
        const flat_index_t begin = offsets[c];
        const flat_index_t end = offsets[c + 1];
        for(flat_index_t k = begin + 1; k < end; ++k)
        {
            const flat_index_t face = faceList[k];
            flat_index_t j = k;
            while(j > begin and faceList[j - 1] > face)
            {
                faceList[j] = faceList[j - 1];
                --j;
            }
            faceList[j] = face;
        }
    });
    Kokkos::fence();
    cellFaceOffsets.Download(hostOut.cellFaceOffsets);
    cellFaces.Download(hostOut.cellFaces);
}

FaceBuildReport KokkosGeometryExecutor::BuildRawFaces(const std::vector<flat_index_t> &order, RawFacesView &raw)
{
    if(InSerialGeometryScope())
    {
        return this->host_.BuildRawFaces(order, raw);
    }
    const std::size_t pointCount = this->mesh_.ownedCount;
    this->order_.Upload(order);
    EnsureSize(this->results_, pointCount, "MadVoro::RawFaceResults");
    EnsureSize(this->faceOffsets_, pointCount + 1, "MadVoro::RawFaceOffsets");
    EnsureSize(this->vertexOffsets_, pointCount + 1, "MadVoro::RawVertexOffsets");
    const auto onHost = [&](void)
    {
        const FaceBuildReport report = BuildRawFacesHost(this->hostMesh_->View(), order.data(), this->hostRaw_);
        raw = HostRawFacesView(*this->hostMesh_, this->hostVertices_->data(), this->hostVertices_->size(), order.data(), this->hostRaw_);
        return report;
    };
    if(pointCount == 0)
    {
        return onHost();
    }
    this->LaunchFaceCount(pointCount);

    // The earliest failure in build order decides; a scratch overflow builds the raw faces on the host.
    const auto results = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), Kokkos::subview(this->results_, std::make_pair(std::size_t(0), pointCount)));
    for(std::size_t i = 0; i < pointCount; ++i)
    {
        if(results(i).error == FaceBuildError::None)
        {
            continue;
        }
        if(results(i).error == FaceBuildError::ScratchOverflow)
        {
            return onHost();
        }
        FaceBuildReport report;
        report.result = results(i);
        report.point = order[i];
        return report;
    }
    flat_index_t faceCount = 0, vertexCount = 0;
    Kokkos::deep_copy(faceCount, Kokkos::subview(this->faceOffsets_, pointCount));
    Kokkos::deep_copy(vertexCount, Kokkos::subview(this->vertexOffsets_, pointCount));
    EnsureSize(this->rawCell_, faceCount, "MadVoro::RawFaceCell");
    EnsureSize(this->rawOther_, faceCount, "MadVoro::RawFaceOther");
    EnsureSize(this->rawVertexOffsets_, faceCount + 1, "MadVoro::RawFaceVertexOffsets");
    EnsureSize(this->rawVertices_, vertexCount, "MadVoro::RawFaceVertices");
    Kokkos::deep_copy(Kokkos::subview(this->rawVertexOffsets_, faceCount), vertexCount);
    this->LaunchFaceFill(pointCount);
    raw = RawFacesView();
    raw.points = this->points_.Data();
    raw.hostPoints = this->hostMesh_->points.data();
    raw.vertices = this->vertices_.Data();
    raw.faceCell = this->rawCell_.data();
    raw.faceOther = this->rawOther_.data();
    raw.faceVertexOffsets = this->rawVertexOffsets_.data();
    raw.faceVertices = this->rawVertices_.data();
    raw.order = order.data();
    raw.cellCount = pointCount;
    raw.pointCount = this->hostMesh_->points.size();
    raw.vertexCount = this->hostVertices_->size();
    raw.faceCount = faceCount;
    raw.faceVertexCount = vertexCount;
    raw.onDevice = true;
    return FaceBuildReport();
}

void KokkosGeometryExecutor::Release(void)
{
    *this = KokkosGeometryExecutor();
}

} // namespace Geometry
} // namespace MadVoro
