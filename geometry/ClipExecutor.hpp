// Runs the cell clipping kernels: on the host, or (KokkosClipExecutor) in Kokkos' default execution space
#ifndef MADVORO_GEOMETRY_CLIP_EXECUTOR_HPP
#define MADVORO_GEOMETRY_CLIP_EXECUTOR_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

#include "FlatMesh.hpp"
#include "HostGeometryExecutor.hpp"
#include "kernels/ClipCell.hpp"
#include "kernels/ClipFaces.hpp"

#if defined(MADVORO_WITH_GPU) and not defined(MADVORO_GEOMETRY_ON_HOST)
    #include "KokkosGeometryExecutor.hpp"
#endif // MADVORO_WITH_GPU && !MADVORO_GEOMETRY_ON_HOST

namespace MadVoro
{
namespace Geometry
{

//! \brief Points sorted into a uniform grid of about two points per grid cell
struct ClipGrid
{
    std::vector<flat_index_t> cellStart;
    std::vector<flat_index_t> points;
    FlatPoint low{0.0, 0.0, 0.0};
    double size = 1;
    int dims[3] = {0, 0, 0};

    ClipGridView View(void) const
    {
        ClipGridView view;
        view.cellStart = this->cellStart.data();
        view.points = this->points.data();
        view.low = this->low;
        view.size = this->size;
        for(int a = 0; a < 3; ++a)
        {
            view.dims[a] = this->dims[a];
        }
        return view;
    }
};

//! \brief The grid of points [first, last) of `points`, but those `skip` marks (no grid cells if there are none)
inline ClipGrid MakeClipGrid(const std::vector<FlatPoint> &points, flat_index_t first, flat_index_t last, const std::vector<std::uint8_t> &skip)
{
    ClipGrid grid;
    std::vector<flat_index_t> indices;
    indices.reserve(last - first);
    for(flat_index_t i = first; i < last; ++i)
    {
        if(skip[i] == 0)
        {
            indices.push_back(i);
        }
    }
    const flat_index_t count = indices.size();
    if(count == 0)
    {
        return grid;
    }
    FlatPoint low = points[indices[0]], high = points[indices[0]];
    for(const flat_index_t i : indices)
    {
        low = FlatPoint{std::min(low.x, points[i].x), std::min(low.y, points[i].y), std::min(low.z, points[i].z)};
        high = FlatPoint{std::max(high.x, points[i].x), std::max(high.y, points[i].y), std::max(high.z, points[i].z)};
    }
    const FlatPoint width = high - low;
    const double extent = std::max({width.x, width.y, width.z});
    const double floor = (extent > 0) ? 1e-6 * extent : 1.0;
    const double volume = std::max(width.x, floor) * std::max(width.y, floor) * std::max(width.z, floor);
    grid.size = std::cbrt(2.0 * volume / static_cast<double>(count));
    if(not (grid.size > 0))
    {
        grid.size = 1.0;
    }
    grid.low = low;
    const double widths[3] = {width.x, width.y, width.z};
    flat_index_t cells = 1;
    for(int a = 0; a < 3; ++a)
    {
        grid.dims[a] = static_cast<int>(std::min(1024.0, std::floor(widths[a] / grid.size) + 1));
        cells *= static_cast<flat_index_t>(grid.dims[a]);
    }
    const auto cellOf = [&grid](const FlatPoint &p)
    {
        const double coordinates[3] = {p.x - grid.low.x, p.y - grid.low.y, p.z - grid.low.z};
        int g[3];
        for(int a = 0; a < 3; ++a)
        {
            g[a] = std::clamp(static_cast<int>(coordinates[a] / grid.size), 0, grid.dims[a] - 1);
        }
        return (static_cast<flat_index_t>(g[0]) * grid.dims[1] + g[1]) * grid.dims[2] + g[2];
    };
    grid.cellStart.assign(cells + 1, 0);
    std::vector<flat_index_t> cellOfPoint(count);
    for(flat_index_t i = 0; i < count; ++i)
    {
        cellOfPoint[i] = cellOf(points[indices[i]]);
        ++grid.cellStart[cellOfPoint[i] + 1];
    }
    for(flat_index_t c = 0; c < cells; ++c)
    {
        grid.cellStart[c + 1] += grid.cellStart[c];
    }
    grid.points.resize(count);
    std::vector<flat_index_t> next(grid.cellStart.begin(), grid.cellStart.end() - 1);
    for(flat_index_t i = 0; i < count; ++i)
    {
        grid.points[next[cellOfPoint[i]]++] = indices[i];
    }
    return grid;
}

//! \brief What failed in a Finish, and in which cell
struct ClipFinishReport
{
    ClipFinishError error = ClipFinishError::None;
    flat_index_t cell = FLAT_INVALID_INDEX;
};

namespace clip_detail
{
    // the cells' arrays on the host
    struct HostClipCells
    {
        std::vector<flat_index_t> planes;
        std::vector<std::uint8_t> planeCount;
        std::vector<ClipTriangle> triangles;
        std::vector<FlatPoint> vertices;
        std::vector<std::uint8_t> triangleCount;
        std::vector<double> radius2, known;
        std::vector<ClipStatus> status;

        void Resize(flat_index_t cellCount)
        {
            this->planes.resize(cellCount * CLIP_MAX_PLANES);
            this->planeCount.resize(cellCount);
            this->triangles.resize(cellCount * CLIP_MAX_VERTICES);
            this->vertices.resize(cellCount * CLIP_MAX_VERTICES);
            this->triangleCount.resize(cellCount);
            this->radius2.resize(cellCount);
            this->known.resize(cellCount);
            this->status.resize(cellCount);
        }

        ClipCellsView View(const std::vector<FlatPoint> &points, flat_index_t cellCount, const FlatPoint &cubeLow, const FlatPoint &cubeHigh, const ClipWalls &walls)
        {
            ClipCellsView view;
            view.points = points.data();
            view.cellCount = cellCount;
            view.cubeLow = cubeLow;
            view.cubeHigh = cubeHigh;
            view.walls = walls;
            view.planes = this->planes.data();
            view.planeCount = this->planeCount.data();
            view.triangles = this->triangles.data();
            view.vertices = this->vertices.data();
            view.triangleCount = this->triangleCount.data();
            view.radius2 = this->radius2.data();
            view.known = this->known.data();
            view.status = this->status.data();
            return view;
        }

        void Release(void)
        {
            *this = HostClipCells();
        }
    };
} // namespace clip_detail

namespace clip_detail
{
    // the raw faces of each cell, collected per cell
    struct ClipRawFaceList
    {
        std::vector<flat_index_t> other, vertexCount, vertices;

        void Add(flat_index_t, flat_index_t neighbor, const flat_index_t *loop, flat_index_t count)
        {
            this->other.push_back(neighbor);
            this->vertexCount.push_back(count);
            this->vertices.insert(this->vertices.end(), loop, loop + count);
        }
    };
} // namespace clip_detail

//! \brief The vertices and raw faces of the cells (cells.firstCell is 0) on the host, the cells' faces in `order`
inline ClipFinishReport RawClipCellsHost(const ClipCellsView &cells, const std::vector<flat_index_t> &order, const ClipMirrorsView &mirrors, std::vector<FlatPoint> &vertices, HostRawFaces &raw)
{
    const flat_index_t cellCount = cells.cellCount;
    std::vector<ClipFinishError> errors(cellCount, ClipFinishError::None);
    HostParallelFor(cellCount, [&](flat_index_t c) { errors[c] = CheckClipCell(cells, c); });
    for(flat_index_t c = 0; c < cellCount; ++c)
    {
        if(errors[c] != ClipFinishError::None)
        {
            return ClipFinishReport{errors[c], c};
        }
    }
    // the vertices
    std::vector<flat_index_t> firstVertex(cellCount + 1, 0);
    for(flat_index_t c = 0; c < cellCount; ++c)
    {
        firstVertex[c + 1] = firstVertex[c] + cells.triangleCount[c];
    }
    vertices.assign(firstVertex[cellCount], FlatPoint{0.0, 0.0, 0.0});
    std::vector<std::uint8_t> exact(firstVertex[cellCount], 0);
    HostParallelFor(cellCount, [&](flat_index_t c)
    {
        for(int t = 0; t < cells.triangleCount[c]; ++t)
        {
            if(not ClipVertexPosition<false>(cells, c, t, vertices[firstVertex[c] + t]))
            {
                exact[firstVertex[c] + t] = 1;
            }
        }
    });
    for(flat_index_t c = 0; c < cellCount; ++c)
    {
        for(int t = 0; t < cells.triangleCount[c]; ++t)
        {
            if(exact[firstVertex[c] + t] != 0)
            {
                ClipVertexPosition<true>(cells, c, t, vertices[firstVertex[c] + t]);
            }
        }
    }
    // the raw faces
    std::vector<clip_detail::ClipRawFaceList> faces(cellCount);
    HostParallelFor(cellCount, [&](flat_index_t c)
    {
        FixedClipFaceScratch<CLIP_MAX_VERTICES> scratch;
        errors[c] = BuildClipCellRawFaces(cells, c, firstVertex[c], mirrors, scratch.View(), faces[c]);
    });
    for(flat_index_t c = 0; c < cellCount; ++c)
    {
        if(errors[c] != ClipFinishError::None)
        {
            return ClipFinishReport{errors[c], c};
        }
    }
    std::vector<flat_index_t> faceOffsets(cellCount + 1, 0), vertexOffsets(cellCount + 1, 0);
    for(flat_index_t i = 0; i < cellCount; ++i)
    {
        const clip_detail::ClipRawFaceList &of = faces[order[i]];
        faceOffsets[i + 1] = faceOffsets[i] + of.other.size();
        vertexOffsets[i + 1] = vertexOffsets[i] + of.vertices.size();
    }
    const flat_index_t faceCount = faceOffsets[cellCount];
    raw.faceCell.resize(faceCount);
    raw.faceOther.resize(faceCount);
    raw.faceVertexOffsets.resize(faceCount + 1);
    raw.faceVertices.resize(vertexOffsets[cellCount]);
    raw.faceVertexOffsets[faceCount] = vertexOffsets[cellCount];
    HostParallelFor(cellCount, [&](flat_index_t i)
    {
        const flat_index_t c = order[i];
        const clip_detail::ClipRawFaceList &of = faces[c];
        flat_index_t vertex = vertexOffsets[i], source = 0;
        for(flat_index_t k = 0; k < of.other.size(); ++k)
        {
            const flat_index_t face = faceOffsets[i] + k;
            raw.faceCell[face] = c;
            raw.faceOther[face] = of.other[k];
            raw.faceVertexOffsets[face] = vertex;
            for(flat_index_t v = 0; v < of.vertexCount[k]; ++v)
            {
                raw.faceVertices[vertex + v] = of.vertices[source + v];
            }
            vertex += of.vertexCount[k];
            source += of.vertexCount[k];
        }
    });
    return ClipFinishReport();
}

//! \brief The view of the clipped cells' raw faces on the host
inline RawFacesView ClipRawFacesView(const std::vector<FlatPoint> &points, flat_index_t cellCount, const std::vector<FlatPoint> &vertices, const std::vector<flat_index_t> &order, const HostRawFaces &raw)
{
    RawFacesView view;
    view.points = points.data();
    view.hostPoints = points.data();
    view.vertices = vertices.data();
    view.faceCell = raw.faceCell.data();
    view.faceOther = raw.faceOther.data();
    view.faceVertexOffsets = raw.faceVertexOffsets.data();
    view.faceVertices = raw.faceVertices.data();
    view.order = order.data();
    view.cellCount = cellCount;
    view.pointCount = points.size();
    view.vertexCount = vertices.size();
    view.faceCount = raw.faceCell.size();
    view.faceVertexCount = raw.faceVertices.size();
    view.onDevice = false;
    return view;
}

//! \brief Calls of a build: Reset, then Update after points are added, then Finish; HostCells reads the cells after Download
class HostClipExecutor
{
public:
    //! \brief Every cell is the cube clipped by the walls; `points` (the cells', then the ghosts') must outlive the calls
    void Reset(const std::vector<FlatPoint> &points, flat_index_t cellCount, const FlatPoint &cubeLow, const FlatPoint &cubeHigh, const ClipWalls &walls)
    {
        this->points_ = &points;
        this->cellCount_ = cellCount;
        this->cubeLow_ = cubeLow;
        this->cubeHigh_ = cubeHigh;
        this->walls_ = walls;
        this->cells_.Resize(cellCount);
        const ClipCellsView cells = this->HostCells();
        HostParallelFor(cellCount, [&cells](flat_index_t c) { StartClipCell<CLIP_EXACT>(cells, c); });
    }

    //! \brief Clips the cells by the points they may not have been clipped by (`added`: the grid of the points from `firstNew` on, or null)
    void Update(const ClipGrid &all, const ClipGrid *added, flat_index_t firstNew)
    {
        const ClipCellsView cells = this->HostCells();
        const ClipGridView allView = all.View();
        const ClipGridView addedView = (added != nullptr) ? added->View() : ClipGridView();
        const ClipGridView *addedPointer = (added != nullptr) ? &addedView : nullptr;
        HostParallelFor(this->cellCount_, [&](flat_index_t c) { UpdateClipCell<CLIP_EXACT>(cells, c, allView, addedPointer, firstNew); });
    }

    //! \brief Nothing to bring back on the host
    void Download(void)
    {
    }

    //! \brief The vertices and raw faces (valid until the next Reset or Finish); mirrorOf[cell * walls + wall]: the cell's mirror ghost across the wall
    ClipFinishReport Finish(const std::vector<flat_index_t> &order, const std::vector<flat_index_t> &mirrorOf, RawFacesView &raw)
    {
        ClipMirrorsView mirrorsView;
        mirrorsView.mirrorOf = mirrorOf.data();
        mirrorsView.wallCount = this->walls_.count;
        const ClipFinishReport report = RawClipCellsHost(this->HostCells(), order, mirrorsView, this->vertices_, this->raw_);
        raw = ClipRawFacesView(*this->points_, this->cellCount_, this->vertices_, order, this->raw_);
        return report;
    }

    //! \brief The cells the last update clipped on the host after the device, and of them those for a vertex (none here)
    std::size_t Resolved(void) const
    {
        return 0;
    }

    std::size_t ResolvedForVertices(void) const
    {
        return 0;
    }

    std::size_t ResolvedWithExpansions(void) const
    {
        return 0;
    }

    //! \brief Each cell's vertex count, and the squared distance from its point to its farthest vertex
    void Radii(std::vector<std::uint8_t> &vertexCounts, std::vector<double> &radius2)
    {
        const ClipCellsView cells = this->HostCells();
        vertexCounts.assign(cells.triangleCount, cells.triangleCount + this->cellCount_);
        radius2.assign(cells.radius2, cells.radius2 + this->cellCount_);
    }

    //! \brief The vertices (not relative to the points) of cells[k], from offsets[k]
    void Vertices(const std::vector<flat_index_t> &cellList, const std::vector<flat_index_t> &offsets, std::vector<FlatPoint> &positions)
    {
        const ClipCellsView cells = this->HostCells();
        positions.resize(offsets.back());
        HostParallelFor(cellList.size(), [&](flat_index_t k)
        {
            const flat_index_t c = cellList[k];
            for(flat_index_t t = 0; t < offsets[k + 1] - offsets[k]; ++t)
            {
                positions[offsets[k] + t] = cells.vertices[c * CLIP_MAX_VERTICES + t] + cells.points[c];
            }
        });
    }

    ClipCellsView HostCells(void)
    {
        return this->cells_.View(*this->points_, this->cellCount_, this->cubeLow_, this->cubeHigh_, this->walls_);
    }

    void Release(void)
    {
        this->cells_.Release();
        this->vertices_ = std::vector<FlatPoint>();
        this->raw_.Release();
    }

private:
    const std::vector<FlatPoint> *points_ = nullptr;
    flat_index_t cellCount_ = 0;
    FlatPoint cubeLow_{0.0, 0.0, 0.0}, cubeHigh_{0.0, 0.0, 0.0};
    ClipWalls walls_;
    clip_detail::HostClipCells cells_;
    std::vector<FlatPoint> vertices_;
    HostRawFaces raw_;
};

#if defined(MADVORO_WITH_GPU) and not defined(MADVORO_GEOMETRY_ON_HOST)
    //! \brief The same calls as HostClipExecutor, with the cells in Kokkos' default memory space (HostCells needs Download)
    class KokkosClipExecutor
    {
    public:
        void Reset(const std::vector<FlatPoint> &points, flat_index_t cellCount, const FlatPoint &cubeLow, const FlatPoint &cubeHigh, const ClipWalls &walls);

        void Update(const ClipGrid &all, const ClipGrid *added, flat_index_t firstNew);

        //! \brief The cells' planes and vertex triangles to the host (not the vertices)
        void Download(void);

        //! \brief The vertices and raw faces in the kernels' memory space (on the host if a face's loop does not fit their scratch), valid until the next Reset or Finish
        ClipFinishReport Finish(const std::vector<flat_index_t> &order, const std::vector<flat_index_t> &mirrorOf, RawFacesView &raw);

        void Radii(std::vector<std::uint8_t> &vertexCounts, std::vector<double> &radius2);

        void Vertices(const std::vector<flat_index_t> &cellList, const std::vector<flat_index_t> &offsets, std::vector<FlatPoint> &positions);

        std::size_t Resolved(void) const
        {
            return this->resolved_;
        }

        //! \brief Of them, those with a vertex floating point could not place
        std::size_t ResolvedForVertices(void) const
        {
            return this->resolvedVertex_;
        }

        //! \brief The cells the last update clipped again in expansions in the kernels
        std::size_t ResolvedWithExpansions(void) const
        {
            return this->resolvedWithExpansions_;
        }

        ClipCellsView HostCells(void);

        void Release(void);

        // Kernel launches (public for CUDA extended lambdas)
        void LaunchInit(void);
        void LaunchUpdate(const ClipGridView &all, const ClipGridView *added, flat_index_t firstNew);
        void LaunchVertices(const flat_index_t *cellList, flat_index_t count, const flat_index_t *offsets, FlatPoint *positions);
        void LaunchScatter(flat_index_t count);
        void LaunchExpansions(const ClipGridView &all, flat_index_t count);
        void LaunchCheck(void);
        void LaunchScatterVertices(flat_index_t count);
        flat_index_t LaunchVertexPositions(void);
        void LaunchFaceCount(const ClipMirrorsView &mirrors);
        void LaunchFaceFill(const ClipMirrorsView &mirrors);

    private:
        // clips again, in expansions in the kernels, the cells floating point could not decide
        void ResolveWithExpansions(const ClipGridView &all);

        // clips exactly, on the host, the cells the kernels could not decide
        void ResolveOnHost(const ClipGrid &all);

        ClipCellsView DeviceCells(void) const;

        const std::vector<FlatPoint> *points_ = nullptr;
        flat_index_t cellCount_ = 0;
        FlatPoint cubeLow_{0.0, 0.0, 0.0}, cubeHigh_{0.0, 0.0, 0.0};
        ClipWalls walls_;
        device_detail::DeviceArray<FlatPoint> devicePoints_;
        device_detail::DeviceArray<flat_index_t> allStart_, allPoints_, addedStart_, addedPoints_;
        Kokkos::View<flat_index_t *> planes_;
        Kokkos::View<std::uint8_t *> planeCount_, triangleCount_;
        Kokkos::View<ClipTriangle *> triangles_;
        Kokkos::View<FlatPoint *> vertices_;
        Kokkos::View<double *> radius2_, known_;
        Kokkos::View<ClipStatus *> status_;
        device_detail::DeviceArray<flat_index_t> vertexCells_, vertexOffsets_;
        device_detail::DeviceArray<FlatPoint> vertexPositions_;
        clip_detail::HostClipCells host_;
        std::size_t resolved_ = 0, resolvedVertex_ = 0, resolvedWithExpansions_ = 0;
        device_detail::DeviceArray<flat_index_t> expansionCells_;
        device_detail::DeviceArray<flat_index_t> resolveCells_, packedPlanes_;
        device_detail::DeviceArray<std::uint8_t> packedPlaneCount_, packedTriangleCount_;
        device_detail::DeviceArray<ClipTriangle> packedTriangles_;
        device_detail::DeviceArray<FlatPoint> packedVertices_;
        device_detail::DeviceArray<double> packedRadius2_, packedKnown_;
        device_detail::DeviceArray<ClipStatus> packedStatus_;
        // the tessellation
        Kokkos::View<ClipFinishError *> errors_;
        Kokkos::View<std::uint8_t *> exact_;
        Kokkos::View<flat_index_t *> faceOffsets_, faceVertexCounts_;
        device_detail::DeviceArray<flat_index_t> order_, mirrorOf_, firstVertex_;
        std::vector<FlatPoint> hostVertices_;
        device_detail::DeviceArray<FlatPoint> tessellationVertices_;
        device_detail::DeviceArray<flat_index_t> exactVertices_;
        device_detail::DeviceArray<FlatPoint> exactPositions_;
        Kokkos::View<flat_index_t *> rawCell_, rawOther_, rawVertexOffsets_, rawVertices_;
        HostRawFaces hostRaw_;
    };
#endif // MADVORO_WITH_GPU && !MADVORO_GEOMETRY_ON_HOST

} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_CLIP_EXECUTOR_HPP
