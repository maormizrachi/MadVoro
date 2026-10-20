#ifndef MADVORO_GEOMETRY_KERNELS_TETRA_VOLUME_HPP
#define MADVORO_GEOMETRY_KERNELS_TETRA_VOLUME_HPP

#include "../FlatMesh.hpp"
#include "../FlatVoronoi.hpp"
#include "../../portability/Math.hpp"

namespace MadVoro
{
namespace Geometry
{

namespace detail
{
    // Shewchuk's static bound for orient3d, as in utils/Predicates3D_internal.hpp.
    constexpr double ORIENT3D_EPSILON = 1.1102230246251565e-016;
    constexpr double ORIENT3D_ERRBOUND_A = (7.0 + 56.0 * ORIENT3D_EPSILON) * ORIENT3D_EPSILON;
} // namespace detail

/*! \brief |orient3d(p0, p1, p2, p3)| / 6, the volume of the tetrahedron, when the
    floating-point determinant is certified by orient3d's static error bound (the
    value orient3d itself returns then). Returns NeedsHost otherwise, where orient3d
    would continue in exact arithmetic; the host then calls orient3d. */
MADVORO_INLINE_FUNCTION GeometryStatus TetraVolume(const FlatPoint &p0, const FlatPoint &p1, const FlatPoint &p2,
                                                   const FlatPoint &p3, double &volume)
{
    const double adx = p0.x - p3.x;
    const double bdx = p1.x - p3.x;
    const double cdx = p2.x - p3.x;
    const double ady = p0.y - p3.y;
    const double bdy = p1.y - p3.y;
    const double cdy = p2.y - p3.y;
    const double adz = p0.z - p3.z;
    const double bdz = p1.z - p3.z;
    const double cdz = p2.z - p3.z;

    const double bdxcdy = bdx * cdy;
    const double cdxbdy = cdx * bdy;
    const double cdxady = cdx * ady;
    const double adxcdy = adx * cdy;
    const double adxbdy = adx * bdy;
    const double bdxady = bdx * ady;

    const double det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);
    const double permanent = (Portable::Abs(bdxcdy) + Portable::Abs(cdxbdy)) * Portable::Abs(adz) +
                             (Portable::Abs(cdxady) + Portable::Abs(adxcdy)) * Portable::Abs(bdz) +
                             (Portable::Abs(adxbdy) + Portable::Abs(bdxady)) * Portable::Abs(cdz);
    const double errbound = detail::ORIENT3D_ERRBOUND_A * permanent;
    if((det > errbound) || (-det > errbound))
    {
        volume = Portable::Abs(det) / 6.0;
        return GeometryStatus::Ok;
    }
    return GeometryStatus::NeedsHost;
}

} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_KERNELS_TETRA_VOLUME_HPP
