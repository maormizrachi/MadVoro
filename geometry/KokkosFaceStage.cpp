#include "FaceStage.hpp"

#include <string>

#include "kernels/CellMoments.hpp"
#include "../portability/KokkosInitializer.hpp"

namespace MadVoro
{
namespace Geometry
{

namespace
{
    using Policy = Kokkos::RangePolicy<Kokkos::IndexType<std::int64_t>>;

    // a face with more vertices sends the stage to the host
    constexpr int FACE_CAPACITY = 64;

    template<typename T>
    void EnsureSize(Kokkos::View<T *> &view, std::size_t size, const char *label)
    {
        if(view.extent(0) < size)
        {
            view = Kokkos::View<T *>(Kokkos::view_alloc(Kokkos::WithoutInitializing, std::string(label)), size);
        }
    }

    template<typename T>
    void CopyToHost(std::vector<T> &host, const T *device, std::size_t size)
    {
        host.resize(size);
        if(size > 0)
        {
            Kokkos::deep_copy(Kokkos::View<T *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(host.data(), size), Kokkos::View<const T *, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(device, size));
        }
    }

    struct KernelScratch
    {
        flat_index_t flip[FACE_CAPACITY];
        flat_index_t backup[FACE_CAPACITY];
        FlatPoint points[FACE_CAPACITY];
        double diffs[FACE_CAPACITY];
        double areas[FACE_CAPACITY];

        MADVORO_INLINE_FUNCTION FaceScratch View()
        {
            FaceScratch view;
            view.flip = this->flip;
            view.backup = this->backup;
            view.points = this->points;
            view.diffs = this->diffs;
            view.areas = this->areas;
            view.capacity = FACE_CAPACITY;
            return view;
        }
    };
} // namespace

void KokkosFaceStage::LaunchFaces(const RawFacesView &raw)
{
    const RawFacesView view = raw;
    const std::uint8_t *looseAreaCut = this->looseAreaCut_.Data();
    const Kokkos::View<RawFaceResult *> results = this->results_;
    const Kokkos::View<flat_index_t *> sizes = this->sizes_, cleaned = this->cleaned_;
    const Kokkos::View<double *> areas = this->areas_;
    const Kokkos::View<FlatPoint *> cms = this->cms_;
    Kokkos::parallel_for("MadVoro::FinishRawFaces", Policy(0, static_cast<std::int64_t>(raw.faceCount)), KOKKOS_LAMBDA(const std::int64_t i)
    {
        const flat_index_t f = static_cast<flat_index_t>(i);
        KernelScratch scratch;
        results(f) = FinishRawFace(view, f, looseAreaCut, scratch.View(), cleaned.data() + view.faceVertexOffsets[f], sizes(f), areas(f), cms(f));
    });
    Kokkos::fence();
}

flat_index_t KokkosFaceStage::LaunchNumbering(flat_index_t faceCount)
{
    const Kokkos::View<RawFaceResult *> results = this->results_;
    const Kokkos::View<flat_index_t *> sizes = this->sizes_, faceIndex = this->faceIndex_, vertexIndex = this->vertexIndex_;
    // each kept face's index and first vertex, in the raw order
    flat_index_t keptCount = 0;
    Kokkos::parallel_scan("MadVoro::RawFaceIndices", Policy(0, static_cast<std::int64_t>(faceCount)), KOKKOS_LAMBDA(const std::int64_t f, flat_index_t &update, const bool final)
    {
        if(final)
        {
            faceIndex(f) = update;
        }
        update += (results(f) == RawFaceResult::Kept) ? 1 : 0;
    }, keptCount);
    Kokkos::parallel_scan("MadVoro::RawFaceVertexIndices", Policy(0, static_cast<std::int64_t>(faceCount)), KOKKOS_LAMBDA(const std::int64_t f, flat_index_t &update, const bool final)
    {
        if(final)
        {
            vertexIndex(f) = update;
        }
        update += (results(f) == RawFaceResult::Kept) ? sizes(f) : 0;
    });
    Kokkos::fence();
    return keptCount;
}

void KokkosFaceStage::LaunchCompaction(const RawFacesView &raw, const FlatVoronoiView &out)
{
    const RawFacesView view = raw;
    const Kokkos::View<RawFaceResult *> results = this->results_;
    const Kokkos::View<flat_index_t *> sizes = this->sizes_, cleaned = this->cleaned_, faceIndex = this->faceIndex_, vertexIndex = this->vertexIndex_;
    const Kokkos::View<double *> areas = this->areas_;
    const Kokkos::View<FlatPoint *> cms = this->cms_;
    Kokkos::parallel_for("MadVoro::CompactRawFaces", Policy(0, static_cast<std::int64_t>(raw.faceCount)), KOKKOS_LAMBDA(const std::int64_t i)
    {
        const flat_index_t f = static_cast<flat_index_t>(i);
        if(results(f) != RawFaceResult::Kept)
        {
            return;
        }
        const flat_index_t face = faceIndex(f);
        out.faceNeighbor0[face] = view.faceCell[f];
        out.faceNeighbor1[face] = face_stage_detail::OutputIndex(view, view.faceOther[f]);
        out.faceArea[face] = areas(f);
        out.faceCM[face] = cms(f);
        out.faceVertexOffsets[face] = vertexIndex(f);
        for(flat_index_t k = 0; k < sizes(f); ++k)
        {
            out.faceVertices[vertexIndex(f) + k] = cleaned(view.faceVertexOffsets[f] + k);
        }
    });
    Kokkos::fence();
}

void KokkosFaceStage::LaunchCellMoments(const RawFacesView &raw, const ConstFlatVoronoiView &view)
{
    const FlatPoint *points = raw.points;
    const flat_index_t *order = this->order_.Data();
    double *volume = this->cellVolume_.Data();
    FlatPoint *cm = this->cellCM_.Data();
    GeometryStatus *status = this->cellStatus_.Data();
    Kokkos::parallel_for("MadVoro::StageCellMoments", Policy(0, static_cast<std::int64_t>(raw.cellCount)), KOKKOS_LAMBDA(const std::int64_t i)
    {
        const flat_index_t c = order[i];
        status[c] = Geometry::CellMoments(c, points[c], view, CertifiedTetraVolume{}, volume[c], cm[c]);
    });
    Kokkos::fence();
}

FaceStageReport KokkosFaceStage::FinishOnHost(const RawFacesView &raw, const std::vector<std::uint8_t> &looseAreaCut, FlatVoronoi &out)
{
    std::vector<FlatPoint> points, vertices;
    std::vector<flat_index_t> faceCell, faceOther, faceVertexOffsets, faceVertices;
    CopyToHost(points, raw.points, raw.pointCount);
    CopyToHost(vertices, raw.vertices, raw.vertexCount);
    CopyToHost(faceCell, raw.faceCell, raw.faceCount);
    CopyToHost(faceOther, raw.faceOther, raw.faceCount);
    CopyToHost(faceVertexOffsets, raw.faceVertexOffsets, raw.faceCount + 1);
    CopyToHost(faceVertices, raw.faceVertices, raw.faceVertexCount);
    RawFacesView host = raw;
    host.points = points.data();
    host.vertices = vertices.data();
    host.faceCell = faceCell.data();
    host.faceOther = faceOther.data();
    host.faceVertexOffsets = faceVertexOffsets.data();
    host.faceVertices = faceVertices.data();
    host.onDevice = false;
    return FinishRawFacesHost(host, looseAreaCut, out);
}

FaceStageReport KokkosFaceStage::Finish(const RawFacesView &raw, const std::vector<std::uint8_t> &looseAreaCut, FlatVoronoi &out)
{
    KokkosInitializer::Instance();
    out = FlatVoronoi();
    const flat_index_t cellCount = raw.cellCount, faceCount = raw.faceCount;
    this->looseAreaCut_.Upload(looseAreaCut);
    EnsureSize(this->results_, faceCount, "MadVoro::FaceStage::results");
    EnsureSize(this->sizes_, faceCount, "MadVoro::FaceStage::sizes");
    EnsureSize(this->cleaned_, raw.faceVertexCount, "MadVoro::FaceStage::cleaned");
    EnsureSize(this->areas_, faceCount, "MadVoro::FaceStage::areas");
    EnsureSize(this->cms_, faceCount, "MadVoro::FaceStage::cms");
    EnsureSize(this->faceIndex_, faceCount + 1, "MadVoro::FaceStage::faceIndex");
    EnsureSize(this->vertexIndex_, faceCount + 1, "MadVoro::FaceStage::vertexIndex");
    this->LaunchFaces(raw);
    std::vector<RawFaceResult> results(faceCount);
    if(faceCount > 0)
    {
        Kokkos::deep_copy(Kokkos::View<RawFaceResult *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(results.data(), faceCount), Kokkos::subview(this->results_, std::make_pair(std::size_t(0), std::size_t(faceCount))));
    }
    for(flat_index_t f = 0; f < faceCount; ++f)
    {
        if(results[f] == RawFaceResult::ScratchOverflow)
        {
            return this->FinishOnHost(raw, looseAreaCut, out);
        }
    }
    std::vector<flat_index_t> faceCell, faceOther;
    for(flat_index_t f = 0; f < faceCount; ++f)
    {
        if(results[f] != RawFaceResult::Kept and results[f] != RawFaceResult::Dropped)
        {
            CopyToHost(faceCell, raw.faceCell, faceCount);
            CopyToHost(faceOther, raw.faceOther, faceCount);
            return FaceStageReport{results[f], faceCell[f], faceOther[f]};
        }
    }
    const flat_index_t keptCount = (faceCount > 0) ? this->LaunchNumbering(faceCount) : 0;
    flat_index_t keptVertices = 0;
    if(faceCount > 0)
    {
        std::vector<flat_index_t> lastSize, lastIndex;
        CopyToHost(lastSize, this->sizes_.data() + (faceCount - 1), 1);
        CopyToHost(lastIndex, this->vertexIndex_.data() + (faceCount - 1), 1);
        keptVertices = lastIndex[0] + ((results[faceCount - 1] == RawFaceResult::Kept) ? lastSize[0] : 0);
    }
    out.ResizeCellsAndTetras(cellCount, raw.vertexCount);
    out.ResizeFaces(keptCount);
    out.faceVertices.assign(keptVertices, 0);
    out.faceVertexOffsets[keptCount] = keptVertices;
    out.cellOrder.assign(raw.order, raw.order + cellCount);
    this->faceNeighbor0_.Output(out.faceNeighbor0);
    this->faceNeighbor1_.Output(out.faceNeighbor1);
    this->faceArea_.Output(out.faceArea);
    this->faceCM_.Output(out.faceCM);
    this->faceVertices_.Output(out.faceVertices);
    this->faceVertexOffsets_.Upload(out.faceVertexOffsets);
    FlatVoronoiView view;
    view.faceVertexOffsets = this->faceVertexOffsets_.Data();
    view.faceVertices = this->faceVertices_.Data();
    view.faceNeighbor0 = this->faceNeighbor0_.Data();
    view.faceNeighbor1 = this->faceNeighbor1_.Data();
    view.faceArea = this->faceArea_.Data();
    view.faceCM = this->faceCM_.Data();
    view.vertices = const_cast<FlatPoint *>(raw.vertices);
    view.cellCount = cellCount;
    view.faceCount = keptCount;
    view.tetraCount = raw.vertexCount;
    this->LaunchCompaction(raw, view);
    device_detail::BuildCellFaceLists(cellCount, keptCount, this->faceNeighbor0_.Data(), this->faceNeighbor1_.Data(), this->cellFaceOffsets_, this->cellFaces_, out);
    // the cells' volumes and centroids (those floating point cannot certify, on the host)
    this->order_.Upload(out.cellOrder);
    this->cellVolume_.Output(out.cellVolume);
    this->cellCM_.Output(out.cellCM);
    this->cellStatus_.Output(out.cellStatus);
    view.cellFaceOffsets = this->cellFaceOffsets_.Data();
    view.cellFaces = this->cellFaces_.Data();
    view.cellVolume = this->cellVolume_.Data();
    view.cellCM = this->cellCM_.Data();
    view.cellStatus = this->cellStatus_.Data();
    this->LaunchCellMoments(raw, view);
    CopyToHost(out.vertices, raw.vertices, raw.vertexCount);
    this->faceNeighbor0_.Download(out.faceNeighbor0);
    this->faceNeighbor1_.Download(out.faceNeighbor1);
    this->faceArea_.Download(out.faceArea);
    this->faceCM_.Download(out.faceCM);
    this->faceVertices_.Download(out.faceVertices);
    this->faceVertexOffsets_.Download(out.faceVertexOffsets);
    this->cellVolume_.Download(out.cellVolume);
    this->cellCM_.Download(out.cellCM);
    this->cellStatus_.Download(out.cellStatus);
    FlatMeshView mesh;
    mesh.points = raw.hostPoints;
    mesh.ownedCount = cellCount;
    ResolveCellMomentsOnHost(mesh, out);
    return FaceStageReport();
}

} // namespace Geometry
} // namespace MadVoro
