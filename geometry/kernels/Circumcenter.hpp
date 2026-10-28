#ifndef MADVORO_GEOMETRY_KERNELS_CIRCUMCENTER_HPP
#define MADVORO_GEOMETRY_KERNELS_CIRCUMCENTER_HPP

#include "../FlatMesh.hpp"
#include "../FlatVoronoi.hpp"
#include "../../portability/Math.hpp"

namespace MadVoro
{
namespace Geometry
{

namespace detail
{
    MADVORO_INLINE_FUNCTION double Det3(double a00, double a01, double a02, double a10, double a11, double a12, double a20, double a21, double a22)
    {
        return a00 * (a11 * a22 - a12 * a21) + a01 * (a12 * a20 - a10 * a22) + a02 * (a10 * a21 - a11 * a20);
    }

    MADVORO_INLINE_FUNCTION double Norm(const FlatPoint &v)
    {
        return Portable::Sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
    }
} // namespace detail

//! \brief Always writes `center`; NeedsHost (and no radius) when too flat or failing the sanity checks.
MADVORO_INLINE_FUNCTION GeometryStatus Circumcenter(const FlatPoint &p0, const FlatPoint &p1, const FlatPoint &p2, const FlatPoint &p3, FlatPoint &center, double &radius)
{
    FlatPoint v2 = p1 - p0;
    FlatPoint v3 = p2 - p0;
    FlatPoint v4 = p3 - p0;

    const double a = detail::Det3(v2.x, v2.y, v2.z, v3.x, v3.y, v3.z, v4.x, v4.y, v4.z);
    if(Portable::Abs(a) < 100 * Portable::DOUBLE_MIN)
    {
        return GeometryStatus::NeedsHost;
    }

    const double s2 = Dot(v2, v2), s3 = Dot(v3, v3), s4 = Dot(v4, v4);
    const double DDx = detail::Det3(s2, v2.y, v2.z, s3, v3.y, v3.z, s4, v4.y, v4.z);
    const double DDy = -detail::Det3(s2, v2.x, v2.z, s3, v3.x, v3.z, s4, v4.x, v4.z);
    const double DDz = detail::Det3(s2, v2.x, v2.y, s3, v3.x, v3.y, s4, v4.x, v4.y);
    center = FlatPoint{DDx / (2 * a), DDy / (2 * a), DDz / (2 * a)} + p0;
    const double Rres = 0.5 * Portable::Sqrt(DDx * DDx + DDy * DDy + DDz * DDz) / Portable::Abs(a);

    // Sanity checks: the four vertices must be equidistant from the center.
    const double Rcheck0 = detail::Norm(p0 - center);
    v2 = v2 + p0;
    v2 = v2 - center;
    const double Rcheck1 = detail::Norm(v2);
    v3 = v3 + p0;
    v3 = v3 - center;
    const double Rcheck2 = detail::Norm(v3);
    v4 = v4 + p0;
    v4 = v4 - center;
    const double Rcheck3 = detail::Norm(v4);
    const double tol = 1 + 1e-6;
    if(((Rcheck0 + Rcheck1 + Rcheck2 + Rcheck3) * tol < (4 * Rcheck0)) or ((Rcheck0 + Rcheck1 + Rcheck2 + Rcheck3) > (tol * 4 * Rcheck0)))
    {
        return GeometryStatus::NeedsHost;
    }
    if(Rcheck0 > tol * Rres or Rcheck0 * tol < Rres)
    {
        return GeometryStatus::NeedsHost;
    }
    const double a_tol = 1e-8;
    if(Portable::Abs(a) < Rres * Rres * Rres * a_tol)
    {
        return GeometryStatus::NeedsHost;
    }
    radius = Rres;
    return GeometryStatus::Ok;
}

MADVORO_INLINE_FUNCTION void CircumcenterOfTetra(flat_index_t t, const FlatMeshView &mesh, const FlatVoronoiView &voronoi)
{
    const FlatTetra &tetra = mesh.tetras[t];
    double radius = 0;
    voronoi.tetraStatus[t] = Circumcenter(mesh.points[tetra.points[0]], mesh.points[tetra.points[1]], mesh.points[tetra.points[2]], mesh.points[tetra.points[3]], voronoi.vertices[t], radius);
    voronoi.tetraRadius[t] = radius;
}

} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_KERNELS_CIRCUMCENTER_HPP
