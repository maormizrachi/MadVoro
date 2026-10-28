// The circumcenter of a tetrahedron the double precision kernel cannot certify, in 50-digit arithmetic (host only)
#ifndef MADVORO_GEOMETRY_HI_PRECISION_CIRCUMCENTER_HPP
#define MADVORO_GEOMETRY_HI_PRECISION_CIRCUMCENTER_HPP

#include <array>
#include <boost/multiprecision/cpp_dec_float.hpp>

#include "FlatMesh.hpp"

namespace MadVoro
{
namespace Geometry
{

//! \brief Writes the circumcenter; returns the circumradius
inline double HiPrecisionCircumcenter(const std::array<FlatPoint, 4> &points, FlatPoint &center)
{
    using Float = boost::multiprecision::cpp_dec_float_50;
    const auto det = [](const std::array<Float, 9> &m)
    {
        return m[0] * (m[4] * m[8] - m[5] * m[7]) + m[1] * (m[5] * m[6] - m[3] * m[8]) + m[2] * (m[3] * m[7] - m[4] * m[6]);
    };
    std::array<std::array<Float, 3>, 4> v;
    for(int k = 0; k < 4; ++k)
    {
        v[k] = {Float(points[k].x), Float(points[k].y), Float(points[k].z)};
    }
    for(int k = 1; k < 4; ++k)
    {
        for(int c = 0; c < 3; ++c)
        {
            v[k][c] -= v[0][c];
        }
    }
    const Float s2 = v[1][0] * v[1][0] + v[1][1] * v[1][1] + v[1][2] * v[1][2];
    const Float s3 = v[2][0] * v[2][0] + v[2][1] * v[2][1] + v[2][2] * v[2][2];
    const Float s4 = v[3][0] * v[3][0] + v[3][1] * v[3][1] + v[3][2] * v[3][2];
    const Float ba = det({v[1][0], v[1][1], v[1][2], v[2][0], v[2][1], v[2][2], v[3][0], v[3][1], v[3][2]});
    const Float bDx = det({s2, v[1][1], v[1][2], s3, v[2][1], v[2][2], s4, v[3][1], v[3][2]});
    const Float bDy = -det({s2, v[1][0], v[1][2], s3, v[2][0], v[2][2], s4, v[3][0], v[3][2]});
    const Float bDz = det({s2, v[1][0], v[1][1], s3, v[2][0], v[2][1], s4, v[3][0], v[3][1]});
    Float temp = (bDx / (2 * ba) + v[0][0]);
    center.x = temp.convert_to<double>();
    temp = (bDy / (2 * ba) + v[0][1]);
    center.y = temp.convert_to<double>();
    temp = (bDz / (2 * ba) + v[0][2]);
    center.z = temp.convert_to<double>();
    temp = (boost::multiprecision::sqrt(bDx * bDx + bDy * bDy + bDz * bDz) / ba);
    return 0.5 * temp.convert_to<double>();
}

} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_HI_PRECISION_CIRCUMCENTER_HPP
