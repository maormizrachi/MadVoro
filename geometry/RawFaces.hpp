// A local builder's faces before the shared stage cleans them: each face once, with its cell, its other point and its vertex loop
#ifndef MADVORO_GEOMETRY_RAW_FACES_HPP
#define MADVORO_GEOMETRY_RAW_FACES_HPP

#include <cstddef>
#include <vector>

#include "FlatMesh.hpp"
#include "kernels/FaceBuild.hpp"

namespace MadVoro
{
namespace Geometry
{

//! \brief The raw faces and what they refer to, in the memory space of the executor that made them (the arrays stay valid until the builder's next build); the faces of each cell are together, the cells in their build order
struct RawFacesView
{
    const FlatPoint *points = nullptr; // the builder's points: the cells', then the gap, then the ghosts
    const FlatPoint *hostPoints = nullptr; // the same on the host
    const FlatPoint *vertices = nullptr;
    const flat_index_t *faceCell = nullptr; // the cell that built the face
    const flat_index_t *faceOther = nullptr;
    const flat_index_t *faceVertexOffsets = nullptr; // faceCount + 1
    const flat_index_t *faceVertices = nullptr; // indices of vertices, in order around the face
    const flat_index_t *order = nullptr; // on the host: the cells in their build order
    flat_index_t cellCount = 0;
    flat_index_t pointCount = 0;
    flat_index_t vertexCount = 0;
    flat_index_t faceCount = 0;
    flat_index_t faceVertexCount = 0;
    // points [gapBegin, gapBegin + gapSize) must not be a kept face's other point (Delaunay's bounding points); in the output, the indices after them shift down by gapSize
    flat_index_t gapBegin = 0;
    flat_index_t gapSize = 0;
    bool onDevice = false;
};

//! \brief Raw faces on the host
struct HostRawFaces
{
    std::vector<flat_index_t> faceCell, faceOther, faceVertexOffsets, faceVertices;

    void Release(void)
    {
        *this = HostRawFaces();
    }
};

//! \brief What failed in the shared stage, and where
struct FaceStageReport
{
    RawFaceResult result = RawFaceResult::Kept;
    flat_index_t cell = FLAT_INVALID_INDEX;
    flat_index_t other = FLAT_INVALID_INDEX;

    bool Failed(void) const
    {
        return this->result != RawFaceResult::Kept;
    }
};

} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_RAW_FACES_HPP
