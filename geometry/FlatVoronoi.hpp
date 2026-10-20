#ifndef MADVORO_GEOMETRY_FLAT_VORONOI_HPP
#define MADVORO_GEOMETRY_FLAT_VORONOI_HPP

#include <cstddef>
#include <cstdint>
#include <type_traits>
#include <vector>

#include "FlatMesh.hpp"
#include "../portability/KokkosTypes.hpp"

namespace MadVoro
{
namespace Geometry
{

/*! \brief Outcome of computing one item (a tetrahedron's circumcenter, or one cell's
    faces) in a portable kernel. NeedsHost marks an item the kernel could not finish,
    e.g. a near-degenerate tetrahedron that needs extended precision; the host then
    recomputes it with the existing code. */
enum class GeometryStatus : std::uint8_t
{
    Ok = 0,
    NeedsHost = 1
};

/*! \brief View of a flat Voronoi tessellation: raw pointers and sizes only, so it can
    be copied into a kernel. The pointers may refer to host or device memory.
    FlatVoronoiView is writable (for the kernels that fill it) and ConstFlatVoronoiView
    is read-only (for its users); a writable view converts to a read-only one.

    Layout, with C cells (the owned points), F faces and T tetrahedra:
    - cell c has faces cellFaces[cellFaceOffsets[c] .. cellFaceOffsets[c + 1])  (CSR)
    - face f has vertices faceVertices[faceVertexOffsets[f] .. faceVertexOffsets[f + 1]),
      in right-hand order with respect to faceNeighbor0[f]  (CSR)
    - a face vertex is a tetrahedron index t, located at vertices[t] (the circumcenter)
    - faceNeighbor0/1[f] are the two mesh points the face separates
    - per face: faceArea, faceCM;  per cell: cellVolume, cellCM, maxRadius, cellStatus;
      per tetrahedron: vertices, tetraRadius, tetraStatus */
template<bool IsConst>
struct BasicFlatVoronoiView
{
    template<typename T>
    using ptr = std::conditional_t<IsConst, const T *, T *>;

    ptr<flat_index_t> cellFaceOffsets = nullptr;
    ptr<flat_index_t> cellFaces = nullptr;
    ptr<flat_index_t> faceVertexOffsets = nullptr;
    ptr<flat_index_t> faceVertices = nullptr;
    ptr<flat_index_t> faceNeighbor0 = nullptr;
    ptr<flat_index_t> faceNeighbor1 = nullptr;
    ptr<double> faceArea = nullptr;
    ptr<FlatPoint> faceCM = nullptr;
    ptr<double> cellVolume = nullptr;
    ptr<FlatPoint> cellCM = nullptr;
    ptr<double> maxRadius = nullptr;
    ptr<GeometryStatus> cellStatus = nullptr;
    ptr<FlatPoint> vertices = nullptr;
    ptr<double> tetraRadius = nullptr;
    ptr<GeometryStatus> tetraStatus = nullptr;
    flat_index_t cellCount = 0;
    flat_index_t faceCount = 0;
    flat_index_t tetraCount = 0;

    BasicFlatVoronoiView() = default;

    //! \brief Read-only view of a writable one.
    template<bool OtherConst, typename = std::enable_if_t<IsConst && !OtherConst>>
    MADVORO_INLINE_FUNCTION BasicFlatVoronoiView(const BasicFlatVoronoiView<OtherConst> &other)
        : cellFaceOffsets(other.cellFaceOffsets), cellFaces(other.cellFaces),
          faceVertexOffsets(other.faceVertexOffsets), faceVertices(other.faceVertices),
          faceNeighbor0(other.faceNeighbor0), faceNeighbor1(other.faceNeighbor1),
          faceArea(other.faceArea), faceCM(other.faceCM), cellVolume(other.cellVolume),
          cellCM(other.cellCM), maxRadius(other.maxRadius), cellStatus(other.cellStatus),
          vertices(other.vertices), tetraRadius(other.tetraRadius), tetraStatus(other.tetraStatus),
          cellCount(other.cellCount), faceCount(other.faceCount), tetraCount(other.tetraCount)
    {
    }

    MADVORO_INLINE_FUNCTION flat_index_t CellFaceCount(flat_index_t cell) const
    {
        return this->cellFaceOffsets[cell + 1] - this->cellFaceOffsets[cell];
    }

    MADVORO_INLINE_FUNCTION flat_index_t CellFace(flat_index_t cell, flat_index_t k) const
    {
        return this->cellFaces[this->cellFaceOffsets[cell] + k];
    }

    MADVORO_INLINE_FUNCTION flat_index_t FaceVertexCount(flat_index_t face) const
    {
        return this->faceVertexOffsets[face + 1] - this->faceVertexOffsets[face];
    }

    MADVORO_INLINE_FUNCTION flat_index_t FaceVertex(flat_index_t face, flat_index_t k) const
    {
        return this->faceVertices[this->faceVertexOffsets[face] + k];
    }

    //! \brief The neighbor of `cell` across `face`.
    MADVORO_INLINE_FUNCTION flat_index_t OtherSide(flat_index_t face, flat_index_t cell) const
    {
        return (this->faceNeighbor0[face] == cell) ? this->faceNeighbor1[face] : this->faceNeighbor0[face];
    }
};

using FlatVoronoiView = BasicFlatVoronoiView<false>;
using ConstFlatVoronoiView = BasicFlatVoronoiView<true>;

static_assert(std::is_trivially_copyable<FlatVoronoiView>::value, "FlatVoronoiView must be trivially copyable");
static_assert(std::is_trivially_copyable<ConstFlatVoronoiView>::value, "ConstFlatVoronoiView must be trivially copyable");

//! \brief Host storage of a flat Voronoi tessellation (layout as in BasicFlatVoronoiView).
struct FlatVoronoi
{
    std::vector<flat_index_t> cellFaceOffsets;
    std::vector<flat_index_t> cellFaces;
    std::vector<flat_index_t> faceVertexOffsets;
    std::vector<flat_index_t> faceVertices;
    std::vector<flat_index_t> faceNeighbor0;
    std::vector<flat_index_t> faceNeighbor1;
    std::vector<double> faceArea;
    std::vector<FlatPoint> faceCM;
    std::vector<double> cellVolume;
    std::vector<FlatPoint> cellCM;
    std::vector<double> maxRadius;
    std::vector<GeometryStatus> cellStatus;
    std::vector<FlatPoint> vertices;
    std::vector<double> tetraRadius;
    std::vector<GeometryStatus> tetraStatus;
    //! \brief The cells in the order they were built (a Hilbert order), for cache-friendly traversal on the host
    std::vector<flat_index_t> cellOrder;

    flat_index_t CellCount() const { return this->cellVolume.size(); }

    flat_index_t FaceCount() const { return this->faceArea.size(); }

    flat_index_t TetraCount() const { return this->vertices.size(); }

    //! \brief Sizes the per-cell and per-tetrahedron arrays, and the cell offsets.
    void ResizeCellsAndTetras(flat_index_t cellCount, flat_index_t tetraCount)
    {
        this->cellFaceOffsets.assign(cellCount + 1, 0);
        this->cellVolume.assign(cellCount, 0.0);
        this->cellCM.assign(cellCount, FlatPoint{0.0, 0.0, 0.0});
        this->maxRadius.assign(cellCount, 0.0);
        this->cellStatus.assign(cellCount, GeometryStatus::Ok);
        this->vertices.assign(tetraCount, FlatPoint{0.0, 0.0, 0.0});
        this->tetraRadius.assign(tetraCount, 0.0);
        this->tetraStatus.assign(tetraCount, GeometryStatus::Ok);
    }

    //! \brief Sizes the per-face arrays and the face offsets.
    void ResizeFaces(flat_index_t faceCount)
    {
        this->faceVertexOffsets.assign(faceCount + 1, 0);
        this->faceNeighbor0.assign(faceCount, FLAT_INVALID_INDEX);
        this->faceNeighbor1.assign(faceCount, FLAT_INVALID_INDEX);
        this->faceArea.assign(faceCount, 0.0);
        this->faceCM.assign(faceCount, FlatPoint{0.0, 0.0, 0.0});
    }

    FlatVoronoiView View()
    {
        FlatVoronoiView view;
        view.cellFaceOffsets = this->cellFaceOffsets.data();
        view.cellFaces = this->cellFaces.data();
        view.faceVertexOffsets = this->faceVertexOffsets.data();
        view.faceVertices = this->faceVertices.data();
        view.faceNeighbor0 = this->faceNeighbor0.data();
        view.faceNeighbor1 = this->faceNeighbor1.data();
        view.faceArea = this->faceArea.data();
        view.faceCM = this->faceCM.data();
        view.cellVolume = this->cellVolume.data();
        view.cellCM = this->cellCM.data();
        view.maxRadius = this->maxRadius.data();
        view.cellStatus = this->cellStatus.data();
        view.vertices = this->vertices.data();
        view.tetraRadius = this->tetraRadius.data();
        view.tetraStatus = this->tetraStatus.data();
        view.cellCount = this->CellCount();
        view.faceCount = this->FaceCount();
        view.tetraCount = this->TetraCount();
        return view;
    }

    ConstFlatVoronoiView View() const
    {
        return const_cast<FlatVoronoi *>(this)->View();
    }
};

/*! \brief Builds the flat layout of an existing tessellation through its public
    accessors (GetPointNo, GetTotalFacesNumber, GetCellFaces, GetPointsInFace,
    GetFaceNeighbors, GetArea, FaceCM, GetVolume, GetCellCM, GetFacePoints), so it
    works for MadVoro::Voronoi3D and for RICH's Tessellation3D alike.
    Face and vertex numbering is kept. Tetrahedron radii and per-cell max radii are
    not part of that interface and are left zero. */
template<typename TessellationT>
FlatVoronoi MakeFlatVoronoi(const TessellationT &tess)
{
    FlatVoronoi flat;
    const flat_index_t cellCount = tess.GetPointNo();
    const flat_index_t faceCount = tess.GetTotalFacesNumber();
    const auto &centers = tess.GetFacePoints();
    flat.ResizeCellsAndTetras(cellCount, centers.size());
    flat.ResizeFaces(faceCount);

    for(flat_index_t t = 0; t < centers.size(); ++t)
    {
        flat.vertices[t] = ToFlatPoint(centers[t]);
    }

    for(flat_index_t c = 0; c < cellCount; ++c)
    {
        flat.cellFaceOffsets[c + 1] = flat.cellFaceOffsets[c] + tess.GetCellFaces(c).size();
        flat.cellVolume[c] = tess.GetVolume(c);
        flat.cellCM[c] = ToFlatPoint(tess.GetCellCM(c));
    }
    flat.cellFaces.reserve(flat.cellFaceOffsets[cellCount]);
    for(flat_index_t c = 0; c < cellCount; ++c)
    {
        const auto &faces = tess.GetCellFaces(c);
        flat.cellFaces.insert(flat.cellFaces.end(), faces.begin(), faces.end());
    }

    for(flat_index_t f = 0; f < faceCount; ++f)
    {
        flat.faceVertexOffsets[f + 1] = flat.faceVertexOffsets[f] + tess.GetPointsInFace(f).size();
        const auto &neighbors = tess.GetFaceNeighbors(f);
        flat.faceNeighbor0[f] = neighbors.first;
        flat.faceNeighbor1[f] = neighbors.second;
        flat.faceArea[f] = tess.GetArea(f);
        flat.faceCM[f] = ToFlatPoint(tess.FaceCM(f));
    }
    flat.faceVertices.reserve(flat.faceVertexOffsets[faceCount]);
    for(flat_index_t f = 0; f < faceCount; ++f)
    {
        const auto &vertices = tess.GetPointsInFace(f);
        flat.faceVertices.insert(flat.faceVertices.end(), vertices.begin(), vertices.end());
    }
    return flat;
}

} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_FLAT_VORONOI_HPP
