// Exact orientation and in-sphere signs in kernels: Shewchuk's expansion arithmetic in fixed-size arrays (an in-sphere determinant that does not fit is left undecided)
#ifndef MADVORO_GEOMETRY_KERNELS_EXPANSIONS_HPP
#define MADVORO_GEOMETRY_KERNELS_EXPANSIONS_HPP

#include "../FlatMesh.hpp"
#include "../../portability/Math.hpp"

namespace MadVoro
{
namespace Geometry
{

// the sign, or EXPANSION_OVERFLOW if the determinant outgrew its array
constexpr int EXPANSION_OVERFLOW = 2;

namespace expansion_detail
{
    // An expansion: components of increasing magnitude that do not overlap, without zeros (but a lone zero)

    // the largest in-sphere determinant (the lift's products are summed into it one component at a time)
    constexpr int DETERMINANT_CAPACITY = 512;
    // a 3 x 3 minor: 3 products of a difference (2 components) and a 2 x 2 minor (16), at most 192 components (with room for the zero it starts from)
    constexpr int MINOR_CAPACITY = 256;

    MADVORO_INLINE_FUNCTION void TwoSum(double a, double b, double &x, double &y)
    {
        x = a + b;
        const double bVirtual = x - a;
        const double aVirtual = x - bVirtual;
        y = (a - aVirtual) + (b - bVirtual);
    }

    // |a| >= |b|
    MADVORO_INLINE_FUNCTION void FastTwoSum(double a, double b, double &x, double &y)
    {
        x = a + b;
        y = b - (x - a);
    }

    MADVORO_INLINE_FUNCTION void TwoDiff(double a, double b, double &x, double &y)
    {
        x = a - b;
        const double bVirtual = a - x;
        const double aVirtual = x + bVirtual;
        y = (a - aVirtual) + (bVirtual - b);
    }

    MADVORO_INLINE_FUNCTION void TwoProduct(double a, double b, double &x, double &y)
    {
        x = a * b;
        y = Portable::Fma(a, b, -x);
    }

    // a - b exactly (1 or 2 components)
    MADVORO_INLINE_FUNCTION int Difference(double a, double b, double *h)
    {
        double x = 0, y = 0;
        TwoDiff(a, b, x, y);
        if(y != 0)
        {
            h[0] = y;
            h[1] = x;
            return 2;
        }
        h[0] = x;
        return 1;
    }

    // h = e + f (Shewchuk's fast_expansion_sum_zeroelim); h holds elen + flen
    MADVORO_INLINE_FUNCTION int Sum(int elen, const double *e, int flen, const double *f, double *h)
    {
        int eindex = 0, findex = 0, hindex = 0;
        double enow = e[0], fnow = f[0];
        double Q = 0, Qnew = 0, hh = 0;
        if((fnow > enow) == (fnow > -enow))
        {
            Q = enow;
            ++eindex;
            enow = (eindex < elen) ? e[eindex] : 0.0;
        }
        else
        {
            Q = fnow;
            ++findex;
            fnow = (findex < flen) ? f[findex] : 0.0;
        }
        if(eindex < elen and findex < flen)
        {
            if((fnow > enow) == (fnow > -enow))
            {
                FastTwoSum(enow, Q, Qnew, hh);
                ++eindex;
                enow = (eindex < elen) ? e[eindex] : 0.0;
            }
            else
            {
                FastTwoSum(fnow, Q, Qnew, hh);
                ++findex;
                fnow = (findex < flen) ? f[findex] : 0.0;
            }
            Q = Qnew;
            if(hh != 0)
            {
                h[hindex++] = hh;
            }
            while(eindex < elen and findex < flen)
            {
                if((fnow > enow) == (fnow > -enow))
                {
                    TwoSum(Q, enow, Qnew, hh);
                    ++eindex;
                    enow = (eindex < elen) ? e[eindex] : 0.0;
                }
                else
                {
                    TwoSum(Q, fnow, Qnew, hh);
                    ++findex;
                    fnow = (findex < flen) ? f[findex] : 0.0;
                }
                Q = Qnew;
                if(hh != 0)
                {
                    h[hindex++] = hh;
                }
            }
        }
        while(eindex < elen)
        {
            TwoSum(Q, enow, Qnew, hh);
            ++eindex;
            enow = (eindex < elen) ? e[eindex] : 0.0;
            Q = Qnew;
            if(hh != 0)
            {
                h[hindex++] = hh;
            }
        }
        while(findex < flen)
        {
            TwoSum(Q, fnow, Qnew, hh);
            ++findex;
            fnow = (findex < flen) ? f[findex] : 0.0;
            Q = Qnew;
            if(hh != 0)
            {
                h[hindex++] = hh;
            }
        }
        if(Q != 0 or hindex == 0)
        {
            h[hindex++] = Q;
        }
        return hindex;
    }

    // h = e * b (Shewchuk's scale_expansion_zeroelim); h holds 2 * elen
    MADVORO_INLINE_FUNCTION int Scale(int elen, const double *e, double b, double *h)
    {
        double Q = 0, hh = 0, product1 = 0, product0 = 0, sum = 0;
        int hindex = 0;
        TwoProduct(e[0], b, Q, hh);
        if(hh != 0)
        {
            h[hindex++] = hh;
        }
        for(int eindex = 1; eindex < elen; ++eindex)
        {
            TwoProduct(e[eindex], b, product1, product0);
            TwoSum(Q, product0, sum, hh);
            if(hh != 0)
            {
                h[hindex++] = hh;
            }
            FastTwoSum(product1, sum, Q, hh);
            if(hh != 0)
            {
                h[hindex++] = hh;
            }
        }
        if(Q != 0 or hindex == 0)
        {
            h[hindex++] = Q;
        }
        return hindex;
    }

    // h = e * f, for e and f of at most 2 components (h holds 8)
    MADVORO_INLINE_FUNCTION int SmallProduct(int elen, const double *e, int flen, const double *f, double *h)
    {
        double first[4], second[4];
        const int firstl = Scale(elen, e, f[0], first);
        if(flen == 1)
        {
            for(int i = 0; i < firstl; ++i)
            {
                h[i] = first[i];
            }
            return firstl;
        }
        const int secondl = Scale(elen, e, f[1], second);
        return Sum(firstl, first, secondl, second, h);
    }

    // h = a * b - c * d, for factors of at most 2 components (h holds 16)
    MADVORO_INLINE_FUNCTION int Minor2(int alen, const double *a, int blen, const double *b, int clen, const double *c, int dlen, const double *d, double *h)
    {
        double ab[8], cd[8];
        const int abl = SmallProduct(alen, a, blen, b, ab);
        const int cdl = SmallProduct(clen, c, dlen, d, cd);
        for(int i = 0; i < cdl; ++i)
        {
            cd[i] = -cd[i];
        }
        return Sum(abl, ab, cdl, cd, h);
    }

    // h += sign * x * m, for x of at most 2 components and m of at most 16 (h and next hold its length + 64)
    MADVORO_INLINE_FUNCTION int AddProduct(int hlen, double *h, int xlen, const double *x, int mlen, const double *m, double sign, double *next)
    {
        double scaled[32];
        for(int i = 0; i < xlen; ++i)
        {
            const int scaledl = Scale(mlen, m, sign * x[i], scaled);
            hlen = Sum(hlen, h, scaledl, scaled, next);
            for(int k = 0; k < hlen; ++k)
            {
                h[k] = next[k];
            }
        }
        return hlen;
    }

    // h = sa * a * ma + sb * b * mb + sc * c * mc (h holds MINOR_CAPACITY)
    MADVORO_INLINE_FUNCTION int Minor3(int al, const double *a, int mal, const double *ma, double sa, int bl, const double *b, int mbl, const double *mb, double sb, int cl, const double *c, int mcl, const double *mc, double sc, double *h)
    {
        double next[MINOR_CAPACITY];
        h[0] = 0;
        int hl = 1;
        hl = AddProduct(hl, h, al, a, mal, ma, sa, next);
        hl = AddProduct(hl, h, bl, b, mbl, mb, sb, next);
        hl = AddProduct(hl, h, cl, c, mcl, mc, sc, next);
        return hl;
    }

    // det += sign * (x^2 + y^2 + z^2) * minor; false if det would outgrow DETERMINANT_CAPACITY
    MADVORO_INLINE_FUNCTION bool AddLiftTimes(int &detl, double *det, double *next, int xl, const double *x, int yl, const double *y, int zl, const double *z, int minorl, const double *minor, double sign)
    {
        const int lengths[3] = {xl, yl, zl};
        const double *coordinates[3] = {x, y, z};
        double scaled[2 * MINOR_CAPACITY];
        for(int axis = 0; axis < 3; ++axis)
        {
            double square[8];
            const int squarel = SmallProduct(lengths[axis], coordinates[axis], lengths[axis], coordinates[axis], square);
            for(int i = 0; i < squarel; ++i)
            {
                const int scaledl = Scale(minorl, minor, sign * square[i], scaled);
                if(detl + scaledl > DETERMINANT_CAPACITY)
                {
                    return false;
                }
                detl = Sum(detl, det, scaledl, scaled, next);
                for(int k = 0; k < detl; ++k)
                {
                    det[k] = next[k];
                }
            }
        }
        return true;
    }

    MADVORO_INLINE_FUNCTION int SignOf(int elen, const double *e)
    {
        const double top = e[elen - 1];
        return (top > 0) ? 1 : ((top < 0) ? -1 : 0);
    }
} // namespace expansion_detail

//! \brief The exact sign of Shewchuk's orient3d of a, b, c and d (its expansions always fit)
MADVORO_INLINE_FUNCTION int ExpansionOrient3dSign(const FlatPoint &a, const FlatPoint &b, const FlatPoint &c, const FlatPoint &d)
{
    using namespace expansion_detail;
    double adx[2], bdx[2], cdx[2], ady[2], bdy[2], cdy[2], adz[2], bdz[2], cdz[2];
    const int adxl = Difference(a.x, d.x, adx), bdxl = Difference(b.x, d.x, bdx), cdxl = Difference(c.x, d.x, cdx);
    const int adyl = Difference(a.y, d.y, ady), bdyl = Difference(b.y, d.y, bdy), cdyl = Difference(c.y, d.y, cdy);
    const int adzl = Difference(a.z, d.z, adz), bdzl = Difference(b.z, d.z, bdz), cdzl = Difference(c.z, d.z, cdz);
    double bc[16], ca[16], ab[16];
    const int bcl = Minor2(bdxl, bdx, cdyl, cdy, cdxl, cdx, bdyl, bdy, bc);
    const int cal = Minor2(cdxl, cdx, adyl, ady, adxl, adx, cdyl, cdy, ca);
    const int abl = Minor2(adxl, adx, bdyl, bdy, bdxl, bdx, adyl, ady, ab);
    // adz * bc + bdz * ca + cdz * ab
    double det[MINOR_CAPACITY];
    const int detl = Minor3(adzl, adz, bcl, bc, 1.0, bdzl, bdz, cal, ca, 1.0, cdzl, cdz, abl, ab, 1.0, det);
    return SignOf(detl, det);
}

//! \brief The exact sign of Shewchuk's insphere of a, b, c and d against e, or EXPANSION_OVERFLOW
MADVORO_INLINE_FUNCTION int ExpansionInsphereSign(const FlatPoint &a, const FlatPoint &b, const FlatPoint &c, const FlatPoint &d, const FlatPoint &e)
{
    using namespace expansion_detail;
    double aex[2], bex[2], cex[2], dex[2], aey[2], bey[2], cey[2], dey[2], aez[2], bez[2], cez[2], dez[2];
    const int aexl = Difference(a.x, e.x, aex), bexl = Difference(b.x, e.x, bex), cexl = Difference(c.x, e.x, cex), dexl = Difference(d.x, e.x, dex);
    const int aeyl = Difference(a.y, e.y, aey), beyl = Difference(b.y, e.y, bey), ceyl = Difference(c.y, e.y, cey), deyl = Difference(d.y, e.y, dey);
    const int aezl = Difference(a.z, e.z, aez), bezl = Difference(b.z, e.z, bez), cezl = Difference(c.z, e.z, cez), dezl = Difference(d.z, e.z, dez);
    // the 2 x 2 minors of the x and y columns
    double ab[16], bc[16], cd[16], da[16], ac[16], bd[16];
    const int abl = Minor2(aexl, aex, beyl, bey, bexl, bex, aeyl, aey, ab);
    const int bcl = Minor2(bexl, bex, ceyl, cey, cexl, cex, beyl, bey, bc);
    const int cdl = Minor2(cexl, cex, deyl, dey, dexl, dex, ceyl, cey, cd);
    const int dal = Minor2(dexl, dex, aeyl, aey, aexl, aex, deyl, dey, da);
    const int acl = Minor2(aexl, aex, ceyl, cey, cexl, cex, aeyl, aey, ac);
    const int bdl = Minor2(bexl, bex, deyl, dey, dexl, dex, beyl, bey, bd);
    // (dlift * abc - clift * dab) + (blift * cda - alift * bcd), with each 3 x 3 minor in turn
    double minor[MINOR_CAPACITY];
    double det[DETERMINANT_CAPACITY], next[DETERMINANT_CAPACITY];
    det[0] = 0;
    int detl = 1;
    int minorl = Minor3(aezl, aez, bcl, bc, 1.0, bezl, bez, acl, ac, -1.0, cezl, cez, abl, ab, 1.0, minor);
    if(not AddLiftTimes(detl, det, next, dexl, dex, deyl, dey, dezl, dez, minorl, minor, 1.0))
    {
        return EXPANSION_OVERFLOW;
    }
    minorl = Minor3(dezl, dez, abl, ab, 1.0, aezl, aez, bdl, bd, 1.0, bezl, bez, dal, da, 1.0, minor);
    if(not AddLiftTimes(detl, det, next, cexl, cex, ceyl, cey, cezl, cez, minorl, minor, -1.0))
    {
        return EXPANSION_OVERFLOW;
    }
    minorl = Minor3(cezl, cez, dal, da, 1.0, dezl, dez, acl, ac, 1.0, aezl, aez, cdl, cd, 1.0, minor);
    if(not AddLiftTimes(detl, det, next, bexl, bex, beyl, bey, bezl, bez, minorl, minor, 1.0))
    {
        return EXPANSION_OVERFLOW;
    }
    minorl = Minor3(bezl, bez, cdl, cd, 1.0, cezl, cez, bdl, bd, -1.0, dezl, dez, bcl, bc, 1.0, minor);
    if(not AddLiftTimes(detl, det, next, aexl, aex, aeyl, aey, aezl, aez, minorl, minor, -1.0))
    {
        return EXPANSION_OVERFLOW;
    }
    return SignOf(detl, det);
}

} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_KERNELS_EXPANSIONS_HPP
