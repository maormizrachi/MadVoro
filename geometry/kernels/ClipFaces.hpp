// The raw tessellation of clipped cells: checks, the vertices' positions and each cell's raw faces, on the host and in kernels
#ifndef MADVORO_GEOMETRY_KERNELS_CLIP_FACES_HPP
#define MADVORO_GEOMETRY_KERNELS_CLIP_FACES_HPP

#include <cstdint>

#include "ClipCell.hpp"
#include "FaceBuild.hpp"

namespace MadVoro
{
namespace Geometry
{

enum class ClipFinishError : std::uint8_t
{
    None = 0,
    CellOverflow, // the cell has more planes or vertices than the capacity
    CellUnresolved, // the cell has a vertex that was not placed or a removal that was not decided
    ReachesCube, // the cell reaches the bounding cube (its ghosts are missing)
    BrokenLoop, // the vertices around a face do not form one loop
    NoMirror, // a face on a wall has no mirror among the ghosts
    ScratchOverflow // a face has more vertices than the scratch (kernels only: the host finishes the build)
};

//! \brief Each cell's mirror ghost across each wall: mirrorOf[cell * wallCount + wall] (FLAT_INVALID_INDEX if it has none)
struct ClipMirrorsView
{
    const flat_index_t *mirrorOf = nullptr;
    int wallCount = 0;
};

//! \brief Arrays for one face's vertex loop at a time
struct ClipFaceScratch
{
    flat_index_t *fan = nullptr;
    int *loopTriangle = nullptr;
    int *loopFrom = nullptr;
    int *loopTo = nullptr;
    int capacity = 0;
};

template<int Capacity>
struct FixedClipFaceScratch
{
    flat_index_t fan[Capacity];
    int loopTriangle[Capacity];
    int loopFrom[Capacity];
    int loopTo[Capacity];

    MADVORO_INLINE_FUNCTION ClipFaceScratch View()
    {
        ClipFaceScratch view;
        view.fan = this->fan;
        view.loopTriangle = this->loopTriangle;
        view.loopFrom = this->loopFrom;
        view.loopTo = this->loopTo;
        view.capacity = Capacity;
        return view;
    }
};

namespace clip_detail
{
    // whether the cell has a vertex on the plane of point q
    MADVORO_INLINE_FUNCTION bool HasPlaneOf(const ClipCellsView &cells, flat_index_t cell, flat_index_t q)
    {
        const flat_index_t local = cell - cells.firstCell;
        const flat_index_t *planes = cells.planes + local * CLIP_MAX_PLANES;
        const ClipTriangle *triangles = cells.triangles + local * CLIP_MAX_VERTICES;
        int slot = -1;
        for(int k = 0; k < cells.planeCount[local]; ++k)
        {
            if(planes[k] == q)
            {
                slot = k;
            }
        }
        if(slot < 0)
        {
            return false;
        }
        for(int t = 0; t < cells.triangleCount[local]; ++t)
        {
            for(int k = 0; k < 3; ++k)
            {
                if(triangles[t].plane[k] == slot)
                {
                    return true;
                }
            }
        }
        return false;
    }

} // namespace clip_detail

//! \brief Whether the cell clipped and stays inside the cube
MADVORO_INLINE_FUNCTION ClipFinishError CheckClipCell(const ClipCellsView &cells, flat_index_t cell)
{
    const flat_index_t local = cell - cells.firstCell;
    if(cells.status[local] == ClipStatus::Overflow)
    {
        return ClipFinishError::CellOverflow;
    }
    if(cells.status[local] != ClipStatus::Ok)
    {
        return ClipFinishError::CellUnresolved;
    }
    const flat_index_t *planes = cells.planes + local * CLIP_MAX_PLANES;
    const ClipTriangle *triangles = cells.triangles + local * CLIP_MAX_VERTICES;
    for(int t = 0; t < cells.triangleCount[local]; ++t)
    {
        for(int k = 0; k < 3; ++k)
        {
            if(planes[triangles[t].plane[k]] >= CLIP_CUBE_PLANE)
            {
                return ClipFinishError::ReachesCube;
            }
        }
    }
    return ClipFinishError::None;
}

/** \brief The position of the cell's vertex t: the circumcenter of the cell's point and the points of the vertex's planes, from the 4 points in increasing index order (mirrors of the first), the same in every cell that has it. False if floating point cannot place it and `Exact` is false */
template<bool Exact>
MADVORO_INLINE_FUNCTION bool ClipVertexPosition(const ClipCellsView &cells, flat_index_t cell, int t, FlatPoint &position)
{
    const flat_index_t local = cell - cells.firstCell;
    const flat_index_t *planes = cells.planes + local * CLIP_MAX_PLANES;
    const ClipTriangle &triangle = cells.triangles[local * CLIP_MAX_VERTICES + t];
    flat_index_t order[4] = {cell, planes[triangle.plane[0]], planes[triangle.plane[1]], planes[triangle.plane[2]]};
    for(int i = 1; i < 4; ++i)
    {
        for(int j = i; j > 0 and order[j] < order[j - 1]; --j)
        {
            const flat_index_t swap = order[j];
            order[j] = order[j - 1];
            order[j - 1] = swap;
        }
    }
    return clip_detail::CanonicalVertex<Exact>(cells, order, position);
}

/** \brief The raw faces cell c builds: its vertices (global vertex firstVertex on) around each of its planes. A face of two cells is built by the lower one, unless only the other has it; a face on a wall is the cell's with its mirror ghost */
MADVORO_HOST_DEVICE_TEMPLATE
template<typename Sink>
MADVORO_INLINE_FUNCTION ClipFinishError BuildClipCellRawFaces(const ClipCellsView &cells, flat_index_t c, flat_index_t firstVertex, const ClipMirrorsView &mirrors, const ClipFaceScratch &s, Sink &sink)
{
    const flat_index_t local = c - cells.firstCell;
    const ClipTriangle *triangles = cells.triangles + local * CLIP_MAX_VERTICES;
    const flat_index_t *planes = cells.planes + local * CLIP_MAX_PLANES;
    const int triangleCount = cells.triangleCount[local];
    for(int slot = 0; slot < cells.planeCount[local]; ++slot)
    {
        const flat_index_t other = planes[slot];
        // the plane's vertices, each with the planes before and after it around the face
        int m = 0;
        for(int t = 0; t < triangleCount; ++t)
        {
            for(int k = 0; k < 3; ++k)
            {
                if(triangles[t].plane[k] == slot)
                {
                    if(m >= s.capacity)
                    {
                        return ClipFinishError::ScratchOverflow;
                    }
                    s.loopTriangle[m] = t;
                    s.loopFrom[m] = triangles[t].plane[(k + 1) % 3];
                    s.loopTo[m] = triangles[t].plane[(k + 2) % 3];
                    ++m;
                }
            }
        }
        if(m < 3)
        {
            continue;
        }
        if(other < cells.firstCell + cells.cellCount and other < c and clip_detail::HasPlaneOf(cells, other, c))
        {
            continue;
        }
        // around the face: the next vertex shares the plane after this one
        int current = 0;
        for(int k = 0; k < m; ++k)
        {
            s.fan[k] = firstVertex + static_cast<flat_index_t>(s.loopTriangle[current]);
            int next = -1;
            for(int j = 0; j < m and next < 0; ++j)
            {
                if(s.loopFrom[j] == s.loopTo[current])
                {
                    next = j;
                }
            }
            if(next < 0 or (next == 0 and k + 1 < m) or (next != 0 and k + 1 == m))
            {
                return ClipFinishError::BrokenLoop;
            }
            current = next;
        }
        const flat_index_t neighbor = (other >= CLIP_WALL_PLANE) ? mirrors.mirrorOf[c * mirrors.wallCount + (other - CLIP_WALL_PLANE)] : other;
        if(neighbor == FLAT_INVALID_INDEX)
        {
            return ClipFinishError::NoMirror;
        }
        sink.Add(c, neighbor, s.fan, static_cast<flat_index_t>(m));
    }
    return ClipFinishError::None;
}

} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_KERNELS_CLIP_FACES_HPP
