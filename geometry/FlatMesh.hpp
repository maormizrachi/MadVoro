#ifndef MADVORO_GEOMETRY_FLAT_MESH_HPP
#define MADVORO_GEOMETRY_FLAT_MESH_HPP

#include <cstddef>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>

#include "../portability/KokkosTypes.hpp"

namespace MadVoro
{
namespace Geometry
{

// Matches the Delaunay's own indices; FLAT_INVALID_INDEX marks removed tetrahedra there too.
using flat_index_t = std::size_t;
constexpr flat_index_t FLAT_INVALID_INDEX = std::numeric_limits<flat_index_t>::max();

struct FlatPoint
{
    double x;
    double y;
    double z;
};

MADVORO_INLINE_FUNCTION FlatPoint operator+(const FlatPoint &a, const FlatPoint &b)
{
    return FlatPoint{a.x + b.x, a.y + b.y, a.z + b.z};
}

MADVORO_INLINE_FUNCTION FlatPoint operator-(const FlatPoint &a, const FlatPoint &b)
{
    return FlatPoint{a.x - b.x, a.y - b.y, a.z - b.z};
}

MADVORO_INLINE_FUNCTION FlatPoint operator*(double s, const FlatPoint &a)
{
    return FlatPoint{s * a.x, s * a.y, s * a.z};
}

MADVORO_INLINE_FUNCTION double Dot(const FlatPoint &a, const FlatPoint &b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

MADVORO_INLINE_FUNCTION FlatPoint Cross(const FlatPoint &a, const FlatPoint &b)
{
    return FlatPoint{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

template<typename PointT>
inline FlatPoint ToFlatPoint(const PointT &point)
{
    return FlatPoint{point.x, point.y, point.z};
}

//! \brief A removed tetrahedron has every entry FLAT_INVALID_INDEX.
struct FlatTetra
{
    flat_index_t points[4];
    flat_index_t neighbors[4];
};

static_assert(std::is_trivially_copyable<FlatPoint>::value, "FlatPoint must be trivially copyable");
static_assert(std::is_trivially_copyable<FlatTetra>::value, "FlatTetra must be trivially copyable");

//! \brief Points [0, ownedCount) are owned; point-to-tetrahedra lists (CSR) exist for those only.
struct FlatMeshView
{
    const FlatPoint *points = nullptr;
    const FlatTetra *tetras = nullptr;
    const flat_index_t *pointTetraOffsets = nullptr;
    const flat_index_t *pointTetras = nullptr;
    flat_index_t pointCount = 0;
    flat_index_t ownedCount = 0;
    flat_index_t tetraCount = 0;

    MADVORO_INLINE_FUNCTION bool IsAlive(flat_index_t tetra) const
    {
        return this->tetras[tetra].points[0] != FLAT_INVALID_INDEX;
    }

    MADVORO_INLINE_FUNCTION bool IsOwned(flat_index_t point) const
    {
        return point < this->ownedCount;
    }

    MADVORO_INLINE_FUNCTION flat_index_t PointTetraCount(flat_index_t point) const
    {
        return this->pointTetraOffsets[point + 1] - this->pointTetraOffsets[point];
    }

    MADVORO_INLINE_FUNCTION flat_index_t PointTetra(flat_index_t point, flat_index_t k) const
    {
        return this->pointTetras[this->pointTetraOffsets[point] + k];
    }
};

struct FlatMesh
{
    std::vector<FlatPoint> points;
    std::vector<FlatTetra> tetras;
    std::vector<flat_index_t> pointTetraOffsets;
    std::vector<flat_index_t> pointTetras;
    flat_index_t ownedCount = 0;

    FlatMeshView View() const
    {
        FlatMeshView view;
        view.points = this->points.data();
        view.tetras = this->tetras.data();
        view.pointTetraOffsets = this->pointTetraOffsets.data();
        view.pointTetras = this->pointTetras.data();
        view.pointCount = this->points.size();
        view.ownedCount = this->ownedCount;
        view.tetraCount = this->tetras.size();
        return view;
    }
};

namespace detail
{
    template<typename DelaunayT>
    void FillPointsAndTetras(FlatMesh &mesh, const DelaunayT &del)
    {
        mesh.ownedCount = del.Norg_;
        mesh.points.resize(del.points_.size());
        for(std::size_t i = 0; i < del.points_.size(); ++i)
        {
            mesh.points[i] = ToFlatPoint(del.points_[i]);
        }
        mesh.tetras.resize(del.tetras_.size());
        for(std::size_t t = 0; t < del.tetras_.size(); ++t)
        {
            for(int k = 0; k < 4; ++k)
            {
                mesh.tetras[t].points[k] = del.tetras_[t].points[k];
                mesh.tetras[t].neighbors[k] = del.tetras_[t].neighbors[k];
            }
        }
        // Free slots keep stale indices until Voronoi3D::SetPointTetras clears them.
        for(const std::size_t &t : del.empty_tetras_)
        {
            for(int k = 0; k < 4; ++k)
            {
                mesh.tetras[t].points[k] = FLAT_INVALID_INDEX;
                mesh.tetras[t].neighbors[k] = FLAT_INVALID_INDEX;
            }
        }
    }
} // namespace detail

//! \brief Builds the point-to-tetrahedra lists in increasing tetrahedron order.
template<typename DelaunayT>
FlatMesh MakeFlatMesh(const DelaunayT &del)
{
    FlatMesh mesh;
    detail::FillPointsAndTetras(mesh, del);
    const flat_index_t ownedCount = mesh.ownedCount;

    mesh.pointTetraOffsets.assign(ownedCount + 1, 0);
    for(const FlatTetra &tetra : mesh.tetras)
    {
        if(tetra.points[0] == FLAT_INVALID_INDEX)
        {
            continue;
        }
        for(int k = 0; k < 4; ++k)
        {
            if(tetra.points[k] < ownedCount)
            {
                ++mesh.pointTetraOffsets[tetra.points[k] + 1];
            }
        }
    }
    for(flat_index_t p = 0; p < ownedCount; ++p)
    {
        mesh.pointTetraOffsets[p + 1] += mesh.pointTetraOffsets[p];
    }

    mesh.pointTetras.resize(mesh.pointTetraOffsets[ownedCount]);
    std::vector<flat_index_t> next(mesh.pointTetraOffsets.begin(), mesh.pointTetraOffsets.end() - 1);
    for(flat_index_t t = 0; t < mesh.tetras.size(); ++t)
    {
        const FlatTetra &tetra = mesh.tetras[t];
        if(tetra.points[0] == FLAT_INVALID_INDEX)
        {
            continue;
        }
        for(int k = 0; k < 4; ++k)
        {
            if(tetra.points[k] < ownedCount)
            {
                mesh.pointTetras[next[tetra.points[k]]++] = t;
            }
        }
    }
    return mesh;
}

//! \brief Keeps the given point-to-tetrahedra lists in their order (it decides the face numbering).
template<typename DelaunayT, typename TetraListT>
FlatMesh MakeFlatMesh(const DelaunayT &del, const std::vector<TetraListT> &pointTetras)
{
    FlatMesh mesh;
    detail::FillPointsAndTetras(mesh, del);
    const flat_index_t ownedCount = mesh.ownedCount;

    mesh.pointTetraOffsets.assign(ownedCount + 1, 0);
    for(flat_index_t p = 0; p < ownedCount; ++p)
    {
        mesh.pointTetraOffsets[p + 1] = mesh.pointTetraOffsets[p] + pointTetras[p].size();
    }
    mesh.pointTetras.reserve(mesh.pointTetraOffsets[ownedCount]);
    for(flat_index_t p = 0; p < ownedCount; ++p)
    {
        mesh.pointTetras.insert(mesh.pointTetras.end(), pointTetras[p].begin(), pointTetras[p].end());
    }
    return mesh;
}

} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_FLAT_MESH_HPP
