// The shared tessellation stage: a local builder's raw faces in, the finished cells out (cleaned faces, their moments, each cell's faces, cell moments), on the host or in Kokkos' default memory space
#ifndef MADVORO_GEOMETRY_FACE_STAGE_HPP
#define MADVORO_GEOMETRY_FACE_STAGE_HPP

#include <algorithm>
#include <cstdint>
#include <vector>

#include "FlatVoronoi.hpp"
#include "HostGeometryExecutor.hpp"
#include "RawFaces.hpp"

#if defined(MADVORO_WITH_GPU) and not defined(MADVORO_GEOMETRY_ON_HOST)
    #include "KokkosGeometryExecutor.hpp"
#endif // MADVORO_WITH_GPU && !MADVORO_GEOMETRY_ON_HOST

namespace MadVoro
{
namespace Geometry
{

namespace face_stage_detail
{
    // a raw face's other point in the output (after the gap)
    MADVORO_INLINE_FUNCTION flat_index_t OutputIndex(const RawFacesView &raw, flat_index_t point)
    {
        return (point >= raw.gapBegin + raw.gapSize) ? point - raw.gapSize : point;
    }

    // scratch for one face at a time on the host
    struct HostScratch
    {
        std::vector<flat_index_t> flip, backup;
        std::vector<FlatPoint> points;
        std::vector<double> diffs, areas;
        FaceScratch view;

        explicit HostScratch(flat_index_t capacity) : flip(capacity), backup(capacity), points(capacity), diffs(capacity), areas(capacity)
        {
            view.flip = this->flip.data();
            view.backup = this->backup.data();
            view.points = this->points.data();
            view.diffs = this->diffs.data();
            view.areas = this->areas.data();
            view.capacity = capacity;
        }
    };

    // each cell's faces, in increasing face order
    inline void CellFaceLists(flat_index_t cellCount, FlatVoronoi &out)
    {
        const flat_index_t faceCount = out.FaceCount();
        std::vector<flat_index_t> &cellOffsets = out.cellFaceOffsets;
        cellOffsets.assign(cellCount + 1, 0);
        for(flat_index_t f = 0; f < faceCount; ++f)
        {
            ++cellOffsets[out.faceNeighbor0[f] + 1];
            if(out.faceNeighbor1[f] < cellCount)
            {
                ++cellOffsets[out.faceNeighbor1[f] + 1];
            }
        }
        for(flat_index_t c = 0; c < cellCount; ++c)
        {
            cellOffsets[c + 1] += cellOffsets[c];
        }
        out.cellFaces.assign(cellOffsets[cellCount], 0);
        std::vector<flat_index_t> next(cellOffsets.begin(), cellOffsets.end() - 1);
        for(flat_index_t f = 0; f < faceCount; ++f)
        {
            out.cellFaces[next[out.faceNeighbor0[f]]++] = f;
            if(out.faceNeighbor1[f] < cellCount)
            {
                out.cellFaces[next[out.faceNeighbor1[f]]++] = f;
            }
        }
    }
} // namespace face_stage_detail

//! \brief The finished cells of raw faces on the host. `looseAreaCut`: per point of the output (no gap), whether its faces take the looser small-area cut
inline FaceStageReport FinishRawFacesHost(const RawFacesView &raw, const std::vector<std::uint8_t> &looseAreaCut, FlatVoronoi &out)
{
    out = FlatVoronoi();
    const flat_index_t cellCount = raw.cellCount, faceCount = raw.faceCount;
    out.ResizeCellsAndTetras(cellCount, raw.vertexCount);
    std::copy(raw.vertices, raw.vertices + raw.vertexCount, out.vertices.begin());
    flat_index_t capacity = 3;
    for(flat_index_t f = 0; f < faceCount; ++f)
    {
        capacity = std::max(capacity, raw.faceVertexOffsets[f + 1] - raw.faceVertexOffsets[f]);
    }
    // each face cleaned into its raw slot, in chunks that each have their scratch
    std::vector<RawFaceResult> results(faceCount);
    std::vector<flat_index_t> sizes(faceCount, 0), cleaned(raw.faceVertexCount);
    std::vector<double> areas(faceCount, 0);
    std::vector<FlatPoint> cms(faceCount);
    const flat_index_t chunks = std::min<flat_index_t>(faceCount, 1024);
    HostParallelFor(chunks, [&](flat_index_t chunk)
    {
        face_stage_detail::HostScratch scratch(capacity);
        for(flat_index_t f = chunk * faceCount / chunks; f < (chunk + 1) * faceCount / chunks; ++f)
        {
            results[f] = FinishRawFace(raw, f, looseAreaCut.data(), scratch.view, cleaned.data() + raw.faceVertexOffsets[f], sizes[f], areas[f], cms[f]);
        }
    });
    for(flat_index_t f = 0; f < faceCount; ++f)
    {
        if(results[f] != RawFaceResult::Kept and results[f] != RawFaceResult::Dropped)
        {
            return FaceStageReport{results[f], raw.faceCell[f], raw.faceOther[f]};
        }
    }
    // the kept faces, numbered in the raw order
    std::vector<flat_index_t> faceIndex(faceCount + 1, 0), vertexIndex(faceCount + 1, 0);
    for(flat_index_t f = 0; f < faceCount; ++f)
    {
        const bool kept = (results[f] == RawFaceResult::Kept);
        faceIndex[f + 1] = faceIndex[f] + (kept ? 1 : 0);
        vertexIndex[f + 1] = vertexIndex[f] + (kept ? sizes[f] : 0);
    }
    const flat_index_t keptCount = faceIndex[faceCount];
    out.ResizeFaces(keptCount);
    out.faceVertices.assign(vertexIndex[faceCount], 0);
    out.faceVertexOffsets[keptCount] = vertexIndex[faceCount];
    out.cellOrder.assign(raw.order, raw.order + cellCount);
    HostParallelFor(faceCount, [&](flat_index_t f)
    {
        if(results[f] != RawFaceResult::Kept)
        {
            return;
        }
        const flat_index_t face = faceIndex[f];
        out.faceNeighbor0[face] = raw.faceCell[f];
        out.faceNeighbor1[face] = face_stage_detail::OutputIndex(raw, raw.faceOther[f]);
        out.faceArea[face] = areas[f];
        out.faceCM[face] = cms[f];
        out.faceVertexOffsets[face] = vertexIndex[f];
        std::copy(cleaned.begin() + raw.faceVertexOffsets[f], cleaned.begin() + raw.faceVertexOffsets[f] + sizes[f], out.faceVertices.begin() + vertexIndex[f]);
    });
    face_stage_detail::CellFaceLists(cellCount, out);
    FlatMeshView mesh;
    mesh.points = raw.points;
    mesh.ownedCount = cellCount;
    CellMomentsHost(mesh, out);
    return FaceStageReport();
}

#if defined(MADVORO_WITH_GPU) and not defined(MADVORO_GEOMETRY_ON_HOST)
    //! \brief The shared stage in Kokkos' default memory space, for raw faces made there (the host finishes a face larger than the kernels' scratch, and the moments floating point cannot certify)
    class KokkosFaceStage
    {
    public:
        FaceStageReport Finish(const RawFacesView &raw, const std::vector<std::uint8_t> &looseAreaCut, FlatVoronoi &out);

        // Kernel launches (public for CUDA extended lambdas)
        void LaunchFaces(const RawFacesView &raw);
        flat_index_t LaunchNumbering(flat_index_t faceCount);
        void LaunchCompaction(const RawFacesView &raw, const FlatVoronoiView &out);
        void LaunchCellMoments(const RawFacesView &raw, const ConstFlatVoronoiView &view);

    private:
        // the raw faces to the host, finished there
        FaceStageReport FinishOnHost(const RawFacesView &raw, const std::vector<std::uint8_t> &looseAreaCut, FlatVoronoi &out);

        device_detail::DeviceArray<std::uint8_t> looseAreaCut_;
        device_detail::DeviceArray<flat_index_t> order_;
        Kokkos::View<RawFaceResult *> results_;
        Kokkos::View<flat_index_t *> sizes_, cleaned_, faceIndex_, vertexIndex_;
        Kokkos::View<double *> areas_;
        Kokkos::View<FlatPoint *> cms_;
        device_detail::DeviceArray<flat_index_t> cellFaceOffsets_, cellFaces_, faceVertexOffsets_, faceVertices_, faceNeighbor0_, faceNeighbor1_;
        device_detail::DeviceArray<double> faceArea_, cellVolume_;
        device_detail::DeviceArray<FlatPoint> faceCM_, cellCM_, vertices_;
        device_detail::DeviceArray<GeometryStatus> cellStatus_;
    };
#endif // MADVORO_WITH_GPU && !MADVORO_GEOMETRY_ON_HOST

//! \brief The shared stage, where the raw faces are
class FaceStage
{
public:
    FaceStageReport Finish(const RawFacesView &raw, const std::vector<std::uint8_t> &looseAreaCut, FlatVoronoi &out)
    {
        #if defined(MADVORO_WITH_GPU) and not defined(MADVORO_GEOMETRY_ON_HOST)
            if(raw.onDevice)
            {
                return this->device_.Finish(raw, looseAreaCut, out);
            }
        #endif // MADVORO_WITH_GPU && !MADVORO_GEOMETRY_ON_HOST
        return FinishRawFacesHost(raw, looseAreaCut, out);
    }

private:
    #if defined(MADVORO_WITH_GPU) and not defined(MADVORO_GEOMETRY_ON_HOST)
        KokkosFaceStage device_;
    #endif // MADVORO_WITH_GPU && !MADVORO_GEOMETRY_ON_HOST
};

} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_FACE_STAGE_HPP
