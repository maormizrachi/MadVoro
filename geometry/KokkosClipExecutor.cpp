#include "ClipExecutor.hpp"

#include <string>
#include <type_traits>

#include "../portability/KokkosInitializer.hpp"

namespace MadVoro
{
namespace Geometry
{

namespace
{
    using Policy = Kokkos::RangePolicy<Kokkos::IndexType<std::int64_t>>;
    // in host memory, the kernels decide exactly themselves
    constexpr bool HOST_ACCESSIBLE = Kokkos::SpaceAccessibility<Kokkos::DefaultExecutionSpace, Kokkos::HostSpace>::accessible;
    constexpr int KERNEL_ARITHMETIC = HOST_ACCESSIBLE ? CLIP_EXACT : CLIP_FLOATING;
    // a face with more vertices sends the tessellation to the host
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
    void CopyToHost(std::vector<T> &host, const Kokkos::View<T *> &device, std::size_t size)
    {
        host.resize(size);
        if(size > 0)
        {
            Kokkos::deep_copy(Kokkos::View<T *, Kokkos::HostSpace, Kokkos::MemoryTraits<Kokkos::Unmanaged>>(host.data(), size), Kokkos::subview(device, std::make_pair(std::size_t(0), size)));
        }
    }
} // namespace

ClipCellsView KokkosClipExecutor::DeviceCells(void) const
{
    ClipCellsView view;
    view.points = this->devicePoints_.Data();
    view.cellCount = this->cellCount_;
    view.cubeLow = this->cubeLow_;
    view.cubeHigh = this->cubeHigh_;
    view.walls = this->walls_;
    view.planes = this->planes_.data();
    view.planeCount = this->planeCount_.data();
    view.triangles = this->triangles_.data();
    view.vertices = this->vertices_.data();
    view.triangleCount = this->triangleCount_.data();
    view.radius2 = this->radius2_.data();
    view.known = this->known_.data();
    view.status = this->status_.data();
    return view;
}

void KokkosClipExecutor::Reset(const std::vector<FlatPoint> &points, flat_index_t cellCount, const FlatPoint &cubeLow, const FlatPoint &cubeHigh, const ClipWalls &walls)
{
    KokkosInitializer::Instance();
    this->points_ = &points;
    this->cellCount_ = cellCount;
    this->cubeLow_ = cubeLow;
    this->cubeHigh_ = cubeHigh;
    this->walls_ = walls;
    EnsureSize(this->planes_, cellCount * CLIP_MAX_PLANES, "MadVoro::Clip::planes");
    EnsureSize(this->planeCount_, cellCount, "MadVoro::Clip::planeCount");
    EnsureSize(this->triangles_, cellCount * CLIP_MAX_VERTICES, "MadVoro::Clip::triangles");
    EnsureSize(this->vertices_, cellCount * CLIP_MAX_VERTICES, "MadVoro::Clip::vertices");
    EnsureSize(this->triangleCount_, cellCount, "MadVoro::Clip::triangleCount");
    EnsureSize(this->radius2_, cellCount, "MadVoro::Clip::radius2");
    EnsureSize(this->known_, cellCount, "MadVoro::Clip::known");
    EnsureSize(this->status_, cellCount, "MadVoro::Clip::status");
    this->devicePoints_.Upload(points);
    this->LaunchInit();
}

void KokkosClipExecutor::LaunchInit(void)
{
    const ClipCellsView cells = this->DeviceCells();
    Kokkos::parallel_for("MadVoro::InitClipCells", Policy(0, static_cast<std::int64_t>(this->cellCount_)), KOKKOS_LAMBDA(const std::int64_t c)
    {
        StartClipCell<KERNEL_ARITHMETIC>(cells, static_cast<flat_index_t>(c));
    });
    Kokkos::fence();
}

void KokkosClipExecutor::Update(const ClipGrid &all, const ClipGrid *added, flat_index_t firstNew)
{
    // the ghosts added since the last update follow the points
    this->devicePoints_.Upload(*this->points_);
    this->allStart_.Upload(all.cellStart);
    this->allPoints_.Upload(all.points);
    ClipGridView allView = all.View();
    allView.cellStart = this->allStart_.Data();
    allView.points = this->allPoints_.Data();
    ClipGridView addedView;
    if(added != nullptr)
    {
        this->addedStart_.Upload(added->cellStart);
        this->addedPoints_.Upload(added->points);
        addedView = added->View();
        addedView.cellStart = this->addedStart_.Data();
        addedView.points = this->addedPoints_.Data();
    }
    this->LaunchUpdate(allView, (added != nullptr) ? &addedView : nullptr, firstNew);
    if constexpr(not HOST_ACCESSIBLE)
    {
        this->ResolveWithExpansions(allView);
        this->ResolveOnHost(all);
    }
}

void KokkosClipExecutor::ResolveWithExpansions(const ClipGridView &all)
{
    std::vector<ClipStatus> status;
    CopyToHost(status, this->status_, this->cellCount_);
    std::vector<flat_index_t> resolve;
    for(flat_index_t c = 0; c < this->cellCount_; ++c)
    {
        if(status[c] == ClipStatus::NeedsHost or status[c] == ClipStatus::NeedsHostVertex)
        {
            resolve.push_back(c);
        }
    }
    this->resolvedWithExpansions_ = resolve.size();
    if(resolve.empty())
    {
        return;
    }
    this->expansionCells_.Upload(resolve);
    this->LaunchExpansions(all, resolve.size());
}

void KokkosClipExecutor::LaunchExpansions(const ClipGridView &all, flat_index_t count)
{
    const ClipCellsView cells = this->DeviceCells();
    const flat_index_t *resolve = this->expansionCells_.Data();
    // each cell again from the start, by every point, deciding the ties floating point cannot in expansions
    Kokkos::parallel_for("MadVoro::ClipCellsWithExpansions", Policy(0, static_cast<std::int64_t>(count)), KOKKOS_LAMBDA(const std::int64_t i)
    {
        StartClipCell<CLIP_EXPANSIONS>(cells, resolve[i]);
        UpdateClipCell<CLIP_EXPANSIONS>(cells, resolve[i], all, nullptr, 0);
    });
    Kokkos::fence();
}

void KokkosClipExecutor::ResolveOnHost(const ClipGrid &all)
{
    std::vector<ClipStatus> status;
    CopyToHost(status, this->status_, this->cellCount_);
    std::vector<flat_index_t> resolve;
    this->resolvedVertex_ = 0;
    for(flat_index_t c = 0; c < this->cellCount_; ++c)
    {
        if(status[c] == ClipStatus::NeedsHost or status[c] == ClipStatus::NeedsHostVertex)
        {
            resolve.push_back(c);
            this->resolvedVertex_ += (status[c] == ClipStatus::NeedsHostVertex) ? 1 : 0;
        }
    }
    this->resolved_ = resolve.size();
    if(resolve.empty())
    {
        return;
    }
    // each cell again from the start, by every point, side by side in packed arrays
    clip_detail::HostClipCells packed;
    packed.Resize(resolve.size());
    const ClipGridView allView = all.View();
    HostParallelFor(resolve.size(), [&](flat_index_t i)
    {
        ClipCellsView cell = packed.View(*this->points_, 1, this->cubeLow_, this->cubeHigh_, this->walls_);
        cell.firstCell = resolve[i];
        cell.planes += i * CLIP_MAX_PLANES;
        cell.planeCount += i;
        cell.triangles += i * CLIP_MAX_VERTICES;
        cell.vertices += i * CLIP_MAX_VERTICES;
        cell.triangleCount += i;
        cell.radius2 += i;
        cell.known += i;
        cell.status += i;
        StartClipCell<CLIP_EXACT>(cell, resolve[i]);
        UpdateClipCell<CLIP_EXACT>(cell, resolve[i], allView, nullptr, 0);
    });
    this->resolveCells_.Upload(resolve);
    this->packedPlanes_.Upload(packed.planes);
    this->packedPlaneCount_.Upload(packed.planeCount);
    this->packedTriangles_.Upload(packed.triangles);
    this->packedVertices_.Upload(packed.vertices);
    this->packedTriangleCount_.Upload(packed.triangleCount);
    this->packedRadius2_.Upload(packed.radius2);
    this->packedKnown_.Upload(packed.known);
    this->packedStatus_.Upload(packed.status);
    this->LaunchScatter(resolve.size());
}

void KokkosClipExecutor::LaunchScatter(flat_index_t count)
{
    const ClipCellsView cells = this->DeviceCells();
    const flat_index_t *resolve = this->resolveCells_.Data();
    const flat_index_t *planes = this->packedPlanes_.Data();
    const std::uint8_t *planeCount = this->packedPlaneCount_.Data();
    const ClipTriangle *triangles = this->packedTriangles_.Data();
    const FlatPoint *vertices = this->packedVertices_.Data();
    const std::uint8_t *triangleCount = this->packedTriangleCount_.Data();
    const double *radius2 = this->packedRadius2_.Data();
    const double *known = this->packedKnown_.Data();
    const ClipStatus *status = this->packedStatus_.Data();
    Kokkos::parallel_for("MadVoro::ScatterClipCells", Policy(0, static_cast<std::int64_t>(count)), KOKKOS_LAMBDA(const std::int64_t j)
    {
        const flat_index_t i = static_cast<flat_index_t>(j);
        const flat_index_t c = resolve[i];
        for(int k = 0; k < CLIP_MAX_PLANES; ++k)
        {
            cells.planes[c * CLIP_MAX_PLANES + k] = planes[i * CLIP_MAX_PLANES + k];
        }
        for(int k = 0; k < CLIP_MAX_VERTICES; ++k)
        {
            cells.triangles[c * CLIP_MAX_VERTICES + k] = triangles[i * CLIP_MAX_VERTICES + k];
            cells.vertices[c * CLIP_MAX_VERTICES + k] = vertices[i * CLIP_MAX_VERTICES + k];
        }
        cells.planeCount[c] = planeCount[i];
        cells.triangleCount[c] = triangleCount[i];
        cells.radius2[c] = radius2[i];
        cells.known[c] = known[i];
        cells.status[c] = status[i];
    });
    Kokkos::fence();
}

void KokkosClipExecutor::LaunchUpdate(const ClipGridView &all, const ClipGridView *added, flat_index_t firstNew)
{
    const ClipCellsView cells = this->DeviceCells();
    const bool hasAdded = (added != nullptr);
    const ClipGridView addedView = hasAdded ? *added : ClipGridView();
    Kokkos::parallel_for("MadVoro::UpdateClipCells", Policy(0, static_cast<std::int64_t>(this->cellCount_)), KOKKOS_LAMBDA(const std::int64_t c)
    {
        UpdateClipCell<KERNEL_ARITHMETIC>(cells, static_cast<flat_index_t>(c), all, hasAdded ? &addedView : nullptr, firstNew);
    });
    Kokkos::fence();
}

void KokkosClipExecutor::Download(void)
{
    const flat_index_t cellCount = this->cellCount_;
    CopyToHost(this->host_.planes, this->planes_, cellCount * CLIP_MAX_PLANES);
    CopyToHost(this->host_.planeCount, this->planeCount_, cellCount);
    CopyToHost(this->host_.triangles, this->triangles_, cellCount * CLIP_MAX_VERTICES);
    CopyToHost(this->host_.triangleCount, this->triangleCount_, cellCount);
    CopyToHost(this->host_.radius2, this->radius2_, cellCount);
    CopyToHost(this->host_.known, this->known_, cellCount);
    CopyToHost(this->host_.status, this->status_, cellCount);
}

ClipCellsView KokkosClipExecutor::HostCells(void)
{
    ClipCellsView view = this->host_.View(*this->points_, this->cellCount_, this->cubeLow_, this->cubeHigh_, this->walls_);
    // the vertices relative to the cells' points stay on the device
    view.vertices = nullptr;
    return view;
}

void KokkosClipExecutor::Radii(std::vector<std::uint8_t> &vertexCounts, std::vector<double> &radius2)
{
    CopyToHost(vertexCounts, this->triangleCount_, this->cellCount_);
    CopyToHost(radius2, this->radius2_, this->cellCount_);
}

void KokkosClipExecutor::Vertices(const std::vector<flat_index_t> &cellList, const std::vector<flat_index_t> &offsets, std::vector<FlatPoint> &positions)
{
    positions.resize(offsets.back());
    if(cellList.empty())
    {
        return;
    }
    this->vertexCells_.Upload(cellList);
    this->vertexOffsets_.Upload(offsets);
    this->vertexPositions_.Output(positions);
    this->LaunchVertices(this->vertexCells_.Data(), cellList.size(), this->vertexOffsets_.Data(), this->vertexPositions_.Data());
    this->vertexPositions_.Download(positions);
}

void KokkosClipExecutor::LaunchVertices(const flat_index_t *cellList, flat_index_t count, const flat_index_t *offsets, FlatPoint *positions)
{
    const ClipCellsView cells = this->DeviceCells();
    Kokkos::parallel_for("MadVoro::ClipVertices", Policy(0, static_cast<std::int64_t>(count)), KOKKOS_LAMBDA(const std::int64_t i)
    {
        const flat_index_t c = cellList[i];
        for(flat_index_t t = 0; t < offsets[i + 1] - offsets[i]; ++t)
        {
            positions[offsets[i] + t] = cells.vertices[c * CLIP_MAX_VERTICES + t] + cells.points[c];
        }
    });
    Kokkos::fence();
}

void KokkosClipExecutor::LaunchCheck(void)
{
    const ClipCellsView cells = this->DeviceCells();
    const Kokkos::View<ClipFinishError *> errors = this->errors_;
    Kokkos::parallel_for("MadVoro::CheckClipCells", Policy(0, static_cast<std::int64_t>(this->cellCount_)), KOKKOS_LAMBDA(const std::int64_t c)
    {
        errors(c) = CheckClipCell(cells, static_cast<flat_index_t>(c));
    });
    Kokkos::fence();
}

flat_index_t KokkosClipExecutor::LaunchVertexPositions(void)
{
    const ClipCellsView cells = this->DeviceCells();
    const flat_index_t *firstVertex = this->firstVertex_.Data();
    FlatPoint *positions = this->tessellationVertices_.Data();
    const Kokkos::View<std::uint8_t *> exact = this->exact_;
    flat_index_t inexact = 0;
    Kokkos::parallel_reduce("MadVoro::ClipVertexPositions", Policy(0, static_cast<std::int64_t>(this->cellCount_)), KOKKOS_LAMBDA(const std::int64_t i, flat_index_t &count)
    {
        const flat_index_t c = static_cast<flat_index_t>(i);
        for(int t = 0; t < cells.triangleCount[c]; ++t)
        {
            const bool placed = ClipVertexPosition<HOST_ACCESSIBLE>(cells, c, t, positions[firstVertex[c] + t]);
            exact(firstVertex[c] + t) = placed ? 0 : 1;
            count += placed ? 0 : 1;
        }
    }, inexact);
    Kokkos::fence();
    return inexact;
}

void KokkosClipExecutor::LaunchScatterVertices(flat_index_t count)
{
    const flat_index_t *indices = this->exactVertices_.Data();
    const FlatPoint *values = this->exactPositions_.Data();
    FlatPoint *positions = this->tessellationVertices_.Data();
    Kokkos::parallel_for("MadVoro::ScatterClipVertices", Policy(0, static_cast<std::int64_t>(count)), KOKKOS_LAMBDA(const std::int64_t i)
    {
        positions[indices[i]] = values[i];
    });
    Kokkos::fence();
}

void KokkosClipExecutor::LaunchFaceCount(const ClipMirrorsView &mirrors)
{
    const ClipCellsView cells = this->DeviceCells();
    const flat_index_t *order = this->order_.Data();
    const flat_index_t *firstVertex = this->firstVertex_.Data();
    const Kokkos::View<ClipFinishError *> errors = this->errors_;
    const Kokkos::View<flat_index_t *> faceOffsets = this->faceOffsets_;
    const Kokkos::View<flat_index_t *> vertexOffsets = this->faceVertexCounts_;
    const std::int64_t count = static_cast<std::int64_t>(this->cellCount_);
    // each cell's raw faces and their vertices, at index i + 1 of its build order ...
    Kokkos::parallel_for("MadVoro::ClipRawFaceCount", Policy(0, count), KOKKOS_LAMBDA(const std::int64_t i)
    {
        const flat_index_t c = order[i];
        FixedClipFaceScratch<FACE_CAPACITY> scratch;
        RawFaceCountSink sink;
        errors(i) = BuildClipCellRawFaces(cells, c, firstVertex[c], mirrors, scratch.View(), sink);
        faceOffsets(i + 1) = sink.faces;
        vertexOffsets(i + 1) = sink.vertices;
        if(i == 0)
        {
            faceOffsets(0) = 0;
            vertexOffsets(0) = 0;
        }
    });
    // ... turned into its first face and vertex
    Kokkos::parallel_scan("MadVoro::ClipRawFaceOffsets", Policy(0, count), KOKKOS_LAMBDA(const std::int64_t i, flat_index_t &update, const bool final)
    {
        update += faceOffsets(i + 1);
        if(final)
        {
            faceOffsets(i + 1) = update;
        }
    });
    Kokkos::parallel_scan("MadVoro::ClipRawVertexOffsets", Policy(0, count), KOKKOS_LAMBDA(const std::int64_t i, flat_index_t &update, const bool final)
    {
        update += vertexOffsets(i + 1);
        if(final)
        {
            vertexOffsets(i + 1) = update;
        }
    });
    Kokkos::fence();
}

void KokkosClipExecutor::LaunchFaceFill(const ClipMirrorsView &mirrors)
{
    const ClipCellsView cells = this->DeviceCells();
    const flat_index_t *order = this->order_.Data();
    const flat_index_t *firstVertex = this->firstVertex_.Data();
    const Kokkos::View<flat_index_t *> faceOffsets = this->faceOffsets_;
    const Kokkos::View<flat_index_t *> vertexOffsets = this->faceVertexCounts_;
    flat_index_t *rawCell = this->rawCell_.data(), *rawOther = this->rawOther_.data(), *rawVertexOffsets = this->rawVertexOffsets_.data(), *rawVertices = this->rawVertices_.data();
    // the same raw faces, written from each cell's first face and vertex
    Kokkos::parallel_for("MadVoro::ClipRawFaceFill", Policy(0, static_cast<std::int64_t>(this->cellCount_)), KOKKOS_LAMBDA(const std::int64_t i)
    {
        const flat_index_t c = order[i];
        FixedClipFaceScratch<FACE_CAPACITY> scratch;
        RawFaceFillSink sink{rawCell, rawOther, rawVertexOffsets, rawVertices, faceOffsets(i), vertexOffsets(i)};
        BuildClipCellRawFaces(cells, c, firstVertex[c], mirrors, scratch.View(), sink);
    });
    Kokkos::fence();
}

ClipFinishReport KokkosClipExecutor::Finish(const std::vector<flat_index_t> &order, const std::vector<flat_index_t> &mirrorOf, RawFacesView &raw)
{
    const flat_index_t cellCount = this->cellCount_;
    ClipMirrorsView mirrorsView;
    mirrorsView.mirrorOf = mirrorOf.data();
    mirrorsView.wallCount = this->walls_.count;
    const auto onHost = [&](void)
    {
        this->Download();
        mirrorsView.mirrorOf = mirrorOf.data();
        const ClipFinishReport report = RawClipCellsHost(this->HostCells(), order, mirrorsView, this->hostVertices_, this->hostRaw_);
        raw = ClipRawFacesView(*this->points_, cellCount, this->hostVertices_, order, this->hostRaw_);
        return report;
    };
    EnsureSize(this->errors_, cellCount, "MadVoro::Clip::errors");
    this->LaunchCheck();
    std::vector<ClipFinishError> errors;
    CopyToHost(errors, this->errors_, cellCount);
    for(flat_index_t c = 0; c < cellCount; ++c)
    {
        if(errors[c] != ClipFinishError::None)
        {
            return ClipFinishReport{errors[c], c};
        }
    }
    // the vertices (those floating point cannot place, on the host)
    std::vector<std::uint8_t> counts;
    CopyToHost(counts, this->triangleCount_, cellCount);
    std::vector<flat_index_t> firstVertex(cellCount + 1, 0);
    for(flat_index_t c = 0; c < cellCount; ++c)
    {
        firstVertex[c + 1] = firstVertex[c] + counts[c];
    }
    const flat_index_t vertexCount = firstVertex[cellCount];
    this->hostVertices_.assign(vertexCount, FlatPoint{0.0, 0.0, 0.0});
    this->firstVertex_.Upload(firstVertex);
    this->tessellationVertices_.Output(this->hostVertices_);
    EnsureSize(this->exact_, vertexCount, "MadVoro::Clip::exact");
    if(this->LaunchVertexPositions() > 0)
    {
        std::vector<std::uint8_t> exact;
        CopyToHost(exact, this->exact_, vertexCount);
        this->Download();
        const ClipCellsView cells = this->HostCells();
        std::vector<flat_index_t> indices;
        std::vector<FlatPoint> positions;
        for(flat_index_t c = 0; c < cellCount; ++c)
        {
            for(int t = 0; t < counts[c]; ++t)
            {
                if(exact[firstVertex[c] + t] != 0)
                {
                    FlatPoint position{0.0, 0.0, 0.0};
                    ClipVertexPosition<true>(cells, c, t, position);
                    indices.push_back(firstVertex[c] + t);
                    positions.push_back(position);
                }
            }
        }
        this->exactVertices_.Upload(indices);
        this->exactPositions_.Upload(positions);
        this->LaunchScatterVertices(indices.size());
    }
    // the raw faces, counted, then written
    this->order_.Upload(order);
    this->mirrorOf_.Upload(mirrorOf);
    mirrorsView.mirrorOf = this->mirrorOf_.Data();
    EnsureSize(this->faceOffsets_, cellCount + 1, "MadVoro::Clip::faceOffsets");
    EnsureSize(this->faceVertexCounts_, cellCount + 1, "MadVoro::Clip::faceVertexOffsets");
    this->LaunchFaceCount(mirrorsView);
    CopyToHost(errors, this->errors_, cellCount);
    for(flat_index_t i = 0; i < cellCount; ++i)
    {
        if(errors[i] == ClipFinishError::ScratchOverflow)
        {
            return onHost();
        }
    }
    for(flat_index_t i = 0; i < cellCount; ++i)
    {
        if(errors[i] != ClipFinishError::None)
        {
            return ClipFinishReport{errors[i], order[i]};
        }
    }
    flat_index_t faceCount = 0, faceVertexCount = 0;
    Kokkos::deep_copy(faceCount, Kokkos::subview(this->faceOffsets_, cellCount));
    Kokkos::deep_copy(faceVertexCount, Kokkos::subview(this->faceVertexCounts_, cellCount));
    EnsureSize(this->rawCell_, faceCount, "MadVoro::Clip::rawCell");
    EnsureSize(this->rawOther_, faceCount, "MadVoro::Clip::rawOther");
    EnsureSize(this->rawVertexOffsets_, faceCount + 1, "MadVoro::Clip::rawVertexOffsets");
    EnsureSize(this->rawVertices_, faceVertexCount, "MadVoro::Clip::rawVertices");
    Kokkos::deep_copy(Kokkos::subview(this->rawVertexOffsets_, faceCount), faceVertexCount);
    this->LaunchFaceFill(mirrorsView);
    raw = RawFacesView();
    raw.points = this->devicePoints_.Data();
    raw.hostPoints = this->points_->data();
    raw.vertices = this->tessellationVertices_.Data();
    raw.faceCell = this->rawCell_.data();
    raw.faceOther = this->rawOther_.data();
    raw.faceVertexOffsets = this->rawVertexOffsets_.data();
    raw.faceVertices = this->rawVertices_.data();
    raw.order = order.data();
    raw.cellCount = cellCount;
    raw.pointCount = this->points_->size();
    raw.vertexCount = vertexCount;
    raw.faceCount = faceCount;
    raw.faceVertexCount = faceVertexCount;
    raw.onDevice = true;
    return ClipFinishReport();
}

void KokkosClipExecutor::Release(void)
{
    *this = KokkosClipExecutor();
}

} // namespace Geometry
} // namespace MadVoro
