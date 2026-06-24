/*
 * New boost multiprecision predicates — commented out in favor of the
 * old Shewchuk-style adaptive-precision implementation (Predicates3D_old.cpp).
 */

#if 0

#include "Predicates3D.hpp"
#include <cmath>
#include <algorithm>
#include <boost/multiprecision/cpp_int.hpp>

namespace {

using int256 = boost::multiprecision::int256_t;
using int384 = boost::multiprecision::number<
    boost::multiprecision::cpp_int_backend<384, 384,
        boost::multiprecision::signed_magnitude,
        boost::multiprecision::checked, void>>;

double const epsilon = 1.1102230246251565e-016;
double const o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon;
double const isperrboundA = (16.0 + 224.0 * epsilon) * epsilon;

// --------------- exact double -> integer conversion ---------------
inline void decompose(double v, long long& mantissa, int& exponent)
{
    if (v == 0.0) { mantissa = 0; exponent = 0; return; }
    int exp;
    double frac = std::frexp(v, &exp);
    mantissa = static_cast<long long>(std::ldexp(frac, 53));
    exponent = exp - 53;
}

template<typename T>
inline T to_exact_int(double v, int e_min)
{
    long long m;
    int e;
    decompose(v, m, e);
    int shift = e - e_min;
    if (shift <= 0) return T(m);
    bool neg = m < 0;
    T abs_val = T(neg ? -m : m) << shift;
    return neg ? -abs_val : abs_val;
}

int find_min_exponent(const Vector3D* pts, std::size_t n)
{
    int e_min = 10000;
    for (std::size_t i = 0; i < n; ++i)
    {
        double coords[3] = {pts[i].x, pts[i].y, pts[i].z};
        for (int c = 0; c < 3; ++c)
        {
            if (coords[c] == 0.0) continue;
            int exp;
            std::frexp(coords[c], &exp);
            int actual_exp = exp - 53;
            if (actual_exp < e_min) e_min = actual_exp;
        }
    }
    if (e_min == 10000) e_min = 0;
    return e_min;
}

// --------------- orient3d exact (int256) ---------------
// 3x3 det with 53-bit mantissa fits in ~170 bits.
double orient3d_exact(std::array<Vector3D, 4> const& points)
{
    int e_min = find_min_exponent(points.data(), 4);

    int256 ax = to_exact_int<int256>(points[0].x, e_min);
    int256 ay = to_exact_int<int256>(points[0].y, e_min);
    int256 az = to_exact_int<int256>(points[0].z, e_min);
    int256 bx = to_exact_int<int256>(points[1].x, e_min);
    int256 by = to_exact_int<int256>(points[1].y, e_min);
    int256 bz = to_exact_int<int256>(points[1].z, e_min);
    int256 cx = to_exact_int<int256>(points[2].x, e_min);
    int256 cy = to_exact_int<int256>(points[2].y, e_min);
    int256 cz = to_exact_int<int256>(points[2].z, e_min);
    int256 dx = to_exact_int<int256>(points[3].x, e_min);
    int256 dy = to_exact_int<int256>(points[3].y, e_min);
    int256 dz = to_exact_int<int256>(points[3].z, e_min);

    int256 Adx = ax - dx, Ady = ay - dy, Adz = az - dz;
    int256 Bdx = bx - dx, Bdy = by - dy, Bdz = bz - dz;
    int256 Cdx = cx - dx, Cdy = cy - dy, Cdz = cz - dz;

    int256 det = Adx * (Bdy * Cdz - Bdz * Cdy)
               + Ady * (Bdz * Cdx - Bdx * Cdz)
               + Adz * (Bdx * Cdy - Bdy * Cdx);

    return static_cast<double>(det);
}

// --------------- insphere exact (int384) ---------------
// 4x4 det with lift column, 53-bit mantissa, det is ~280 bits.
double insphere_exact(std::array<Vector3D, 5> const& points)
{
    int e_min = find_min_exponent(points.data(), 5);

    int384 ax = to_exact_int<int384>(points[0].x, e_min);
    int384 ay = to_exact_int<int384>(points[0].y, e_min);
    int384 az = to_exact_int<int384>(points[0].z, e_min);
    int384 bx = to_exact_int<int384>(points[1].x, e_min);
    int384 by = to_exact_int<int384>(points[1].y, e_min);
    int384 bz = to_exact_int<int384>(points[1].z, e_min);
    int384 cx = to_exact_int<int384>(points[2].x, e_min);
    int384 cy = to_exact_int<int384>(points[2].y, e_min);
    int384 cz = to_exact_int<int384>(points[2].z, e_min);
    int384 dxv = to_exact_int<int384>(points[3].x, e_min);
    int384 dyv = to_exact_int<int384>(points[3].y, e_min);
    int384 dzv = to_exact_int<int384>(points[3].z, e_min);
    int384 ex = to_exact_int<int384>(points[4].x, e_min);
    int384 ey = to_exact_int<int384>(points[4].y, e_min);
    int384 ez = to_exact_int<int384>(points[4].z, e_min);

    int384 aex = ax - ex, aey = ay - ey, aez = az - ez;
    int384 bex = bx - ex, bey = by - ey, bez = bz - ez;
    int384 cex = cx - ex, cey = cy - ey, cez = cz - ez;
    int384 dex = dxv - ex, dey = dyv - ey, dez = dzv - ez;

    int384 ab = aex * bey - aey * bex;
    int384 bc = bex * cey - bey * cex;
    int384 cd = cex * dey - cey * dex;
    int384 da = dex * aey - dey * aex;
    int384 ac = aex * cey - aey * cex;
    int384 bd = bex * dey - bey * dex;

    int384 abc = aez * bc - bez * ac + cez * ab;
    int384 bcd = bez * cd - cez * bd + dez * bc;
    int384 cda = cez * da + dez * ac + aez * cd;
    int384 dab = dez * ab + aez * bd + bez * da;

    int384 alift = aex * aex + aey * aey + aez * aez;
    int384 blift = bex * bex + bey * bey + bez * bez;
    int384 clift = cex * cex + cey * cey + cez * cez;
    int384 dlift = dex * dex + dey * dey + dez * dez;

    int384 det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

    return static_cast<double>(det);
}

} // anonymous namespace

// --------------- public API ---------------

double orient3d(std::array<Vector3D, 4> const& points)
{
    double adx = points[0].x - points[3].x;
    double bdx = points[1].x - points[3].x;
    double cdx = points[2].x - points[3].x;
    double ady = points[0].y - points[3].y;
    double bdy = points[1].y - points[3].y;
    double cdy = points[2].y - points[3].y;
    double adz = points[0].z - points[3].z;
    double bdz = points[1].z - points[3].z;
    double cdz = points[2].z - points[3].z;

    double bdxcdy = bdx * cdy;
    double cdxbdy = cdx * bdy;
    double cdxady = cdx * ady;
    double adxcdy = adx * cdy;
    double adxbdy = adx * bdy;
    double bdxady = bdx * ady;

    double det = adz * (bdxcdy - cdxbdy)
               + bdz * (cdxady - adxcdy)
               + cdz * (adxbdy - bdxady);

    double permanent = (std::abs(bdxcdy) + std::abs(cdxbdy)) * std::abs(adz)
                     + (std::abs(cdxady) + std::abs(adxcdy)) * std::abs(bdz)
                     + (std::abs(adxbdy) + std::abs(bdxady)) * std::abs(cdz);
    double errbound = o3derrboundA * permanent;
    if ((det > errbound) || (-det > errbound))
        return det;

    return orient3d_exact(points);
}

double insphere(std::array<Vector3D, 5> const& points)
{
    double aex = points[0].x - points[4].x;
    double bex = points[1].x - points[4].x;
    double cex = points[2].x - points[4].x;
    double dex = points[3].x - points[4].x;
    double aey = points[0].y - points[4].y;
    double bey = points[1].y - points[4].y;
    double cey = points[2].y - points[4].y;
    double dey = points[3].y - points[4].y;
    double aez = points[0].z - points[4].z;
    double bez = points[1].z - points[4].z;
    double cez = points[2].z - points[4].z;
    double dez = points[3].z - points[4].z;

    double aexbey = aex * bey, bexaey = bex * aey;
    double bexcey = bex * cey, cexbey = cex * bey;
    double cexdey = cex * dey, dexcey = dex * cey;
    double dexaey = dex * aey, aexdey = aex * dey;
    double aexcey = aex * cey, cexaey = cex * aey;
    double bexdey = bex * dey, dexbey = dex * bey;

    double ab = aexbey - bexaey;
    double bc = bexcey - cexbey;
    double cd = cexdey - dexcey;
    double da = dexaey - aexdey;
    double ac = aexcey - cexaey;
    double bd = bexdey - dexbey;

    double abc = aez * bc - bez * ac + cez * ab;
    double bcd = bez * cd - cez * bd + dez * bc;
    double cda = cez * da + dez * ac + aez * cd;
    double dab = dez * ab + aez * bd + bez * da;

    double alift = aex * aex + aey * aey + aez * aez;
    double blift = bex * bex + bey * bey + bez * bez;
    double clift = cex * cex + cey * cey + cez * cez;
    double dlift = dex * dex + dey * dey + dez * dez;

    double det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

    double permanent = ((std::abs(cexdey) + std::abs(dexcey)) * std::abs(bez)
                      + (std::abs(dexbey) + std::abs(bexdey)) * std::abs(cez)
                      + (std::abs(bexcey) + std::abs(cexbey)) * std::abs(dez)) * alift
                     + ((std::abs(dexaey) + std::abs(aexdey)) * std::abs(cez)
                      + (std::abs(aexcey) + std::abs(cexaey)) * std::abs(dez)
                      + (std::abs(cexdey) + std::abs(dexcey)) * std::abs(aez)) * blift
                     + ((std::abs(aexbey) + std::abs(bexaey)) * std::abs(dez)
                      + (std::abs(bexdey) + std::abs(dexbey)) * std::abs(aez)
                      + (std::abs(dexaey) + std::abs(aexdey)) * std::abs(bez)) * clift
                     + ((std::abs(bexcey) + std::abs(cexbey)) * std::abs(aez)
                      + (std::abs(cexaey) + std::abs(aexcey)) * std::abs(bez)
                      + (std::abs(aexbey) + std::abs(bexaey)) * std::abs(cez)) * dlift;
    double errbound = isperrboundA * permanent;
    if ((det > errbound) || (-det > errbound))
        return det;

    return insphere_exact(points);
}

#endif
