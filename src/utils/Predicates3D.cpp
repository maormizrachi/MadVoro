#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Predicates3D.hpp"

double const epsilon = 1.1102230246251565e-016;
double const splitter = 134217729;

namespace MadVoro
{
	namespace
	{
	#define Absolute(a)  std::abs(a)//((a) >= 0.0 ? (a) : -(a))

	#define Fast_Two_Sum_Tail(a, b, x, y) \
	bvirt = x - a; \
	y = b - bvirt

	#define Fast_Two_Sum(a, b, x, y) \
	x = static_cast<double>(a + b); \
	Fast_Two_Sum_Tail(a, b, x, y)

	/*
	#define Fast_Two_Diff_Tail(a, b, x, y) \
	bvirt = a - x; \
	y = bvirt - b
	*/

	/*
	#define Fast_Two_Diff(a, b, x, y) \
	x = (double) (a - b); \
	Fast_Two_Diff_Tail(a, b, x, y)
	*/

	#define Two_Sum_Tail(a, b, x, y) \
	bvirt = static_cast<double>(x - a); \
	avirt = x - bvirt; \
	bround = b - bvirt; \
	around = a - avirt; \
	y = around + bround

	#define Two_Sum(a, b, x, y) \
	x = static_cast<double>(a + b); \
	Two_Sum_Tail(a, b, x, y)

	#define Two_Diff_Tail(a, b, x, y) \
	bvirt = static_cast<double>(a - x); \
	avirt = x + bvirt; \
	bround = bvirt - b; \
	around = a - avirt; \
	y = around + bround

	#define Two_Diff(a, b, x, y) \
	x = static_cast<double> (a - b); \
	Two_Diff_Tail(a, b, x, y)

	#define Split(a, ahi, alo) \
	c = static_cast<double>(splitter * a); \
	abig = static_cast<double>(c - a); \
	ahi = c - abig; \
	alo = a - ahi

	#define Two_Product_Tail(a, b, x, y) \
	Split(a, ahi, alo); \
	Split(b, bhi, blo); \
	err1 = x - (ahi * bhi); \
	err2 = err1 - (alo * bhi); \
	err3 = err2 - (ahi * blo); \
	y = (alo * blo) - err3

	#define Two_Product(a, b, x, y) \
	x = static_cast<double> (a * b); \
	Two_Product_Tail(a, b, x, y)

	#define Two_Product_Presplit(a, b, bhi, blo, x, y) \
	x = static_cast<double> (a * b); \
	Split(a, ahi, alo); \
	err1 = x - (ahi * bhi); \
	err2 = err1 - (alo * bhi); \
	err3 = err2 - (ahi * blo); \
	y = (alo * blo) - err3

	/*
	#define Two_Product_2Presplit(a, ahi, alo, b, bhi, blo, x, y) \
	x = (double) (a * b); \
	err1 = x - (ahi * bhi); \
	err2 = err1 - (alo * bhi); \
	err3 = err2 - (ahi * blo); \
	y = (alo * blo) - err3
	*/

	/*
	#define Square_Tail(a, x, y) \
	Split(a, ahi, alo); \
	err1 = x - (ahi * ahi); \
	err3 = err1 - ((ahi + ahi) * alo); \
	y = (alo * alo) - err3
	*/

	/*
	#define Square(a, x, y) \
	x = (double) (a * a); \
	Square_Tail(a, x, y)
	*/

	/*
	#define Two_One_Sum(a1, a0, b, x2, x1, x0) \
	Two_Sum(a0, b , _i, x0); \
	Two_Sum(a1, _i, x2, x1)
	*/

	#define Two_One_Diff(a1, a0, b, x2, x1, x0) \
	Two_Diff(a0, b , _i, x0); \
	Two_Sum( a1, _i, x2, x1)

	/*
	#define Two_Two_Sum(a1, a0, b1, b0, x3, x2, x1, x0) \
	Two_One_Sum(a1, a0, b0, _j, _0, x0); \
	Two_One_Sum(_j, _0, b1, x3, x2, x1)
	*/

	#define Two_Two_Diff(a1, a0, b1, b0, x3, x2, x1, x0) \
	Two_One_Diff(a1, a0, b0, _j, _0, x0); \
	Two_One_Diff(_j, _0, b1, x3, x2, x1)

	/*
	#define Four_One_Sum(a3, a2, a1, a0, b, x4, x3, x2, x1, x0) \
	Two_One_Sum(a1, a0, b , _j, x1, x0); \
	Two_One_Sum(a3, a2, _j, x4, x3, x2)
	*/

	/*
	#define Four_Two_Sum(a3, a2, a1, a0, b1, b0, x5, x4, x3, x2, x1, x0) \
	Four_One_Sum(a3, a2, a1, a0, b0, _k, _2, _1, _0, x0); \
	Four_One_Sum(_k, _2, _1, _0, b1, x5, x4, x3, x2, x1)
	*/

	/*
	#define Four_Four_Sum(a3, a2, a1, a0, b4, b3, b1, b0, x7, x6, x5, x4, x3, x2, \
						x1, x0) \
	Four_Two_Sum(a3, a2, a1, a0, b1, b0, _l, _2, _1, _0, x1, x0); \
	Four_Two_Sum(_l, _2, _1, _0, b4, b3, x7, x6, x5, x4, x3, x2)
	*/

	/*
	#define Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b, x8, x7, x6, x5, x4, \
						x3, x2, x1, x0) \
	Four_One_Sum(a3, a2, a1, a0, b , _j, x3, x2, x1, x0); \
	Four_One_Sum(a7, a6, a5, a4, _j, x8, x7, x6, x5, x4)
	*/

	/*
	#define Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0, x9, x8, x7, \
						x6, x5, x4, x3, x2, x1, x0) \
	Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b0, _k, _6, _5, _4, _3, _2, \
					_1, _0, x0); \
	Eight_One_Sum(_k, _6, _5, _4, _3, _2, _1, _0, b1, x9, x8, x7, x6, x5, x4, \
					x3, x2, x1)
	*/

	/*
	#define Eight_Four_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b4, b3, b1, b0, x11, \
						x10, x9, x8, x7, x6, x5, x4, x3, x2, x1, x0) \
	Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0, _l, _6, _5, _4, _3, \
					_2, _1, _0, x1, x0); \
	Eight_Two_Sum(_l, _6, _5, _4, _3, _2, _1, _0, b4, b3, x11, x10, x9, x8, \
					x7, x6, x5, x4, x3, x2)
	*/

	#define Two_One_Product(a1, a0, b, x3, x2, x1, x0) \
	Split(b, bhi, blo); \
	Two_Product_Presplit(a0, b, bhi, blo, _i, x0); \
	Two_Product_Presplit(a1, b, bhi, blo, _j, _0); \
	Two_Sum(_i, _0, _k, x1); \
	Fast_Two_Sum(_j, _k, x3, x2)

	/*
	#define Four_One_Product(a3, a2, a1, a0, b, x7, x6, x5, x4, x3, x2, x1, x0) \
	Split(b, bhi, blo); \
	Two_Product_Presplit(a0, b, bhi, blo, _i, x0); \
	Two_Product_Presplit(a1, b, bhi, blo, _j, _0); \
	Two_Sum(_i, _0, _k, x1); \
	Fast_Two_Sum(_j, _k, _i, x2); \
	Two_Product_Presplit(a2, b, bhi, blo, _j, _0); \
	Two_Sum(_i, _0, _k, x3); \
	Fast_Two_Sum(_j, _k, _i, x4); \
	Two_Product_Presplit(a3, b, bhi, blo, _j, _0); \
	Two_Sum(_i, _0, _k, x5); \
	Fast_Two_Sum(_j, _k, x7, x6)
	*/

	/*
	#define Two_Two_Product(a1, a0, b1, b0, x7, x6, x5, x4, x3, x2, x1, x0) \
	Split(a0, a0hi, a0lo); \
	Split(b0, bhi, blo); \
	Two_Product_2Presplit(a0, a0hi, a0lo, b0, bhi, blo, _i, x0); \
	Split(a1, a1hi, a1lo); \
	Two_Product_2Presplit(a1, a1hi, a1lo, b0, bhi, blo, _j, _0); \
	Two_Sum(_i, _0, _k, _1); \
	Fast_Two_Sum(_j, _k, _l, _2); \
	Split(b1, bhi, blo); \
	Two_Product_2Presplit(a0, a0hi, a0lo, b1, bhi, blo, _i, _0); \
	Two_Sum(_1, _0, _k, x1); \
	Two_Sum(_2, _k, _j, _1); \
	Two_Sum(_l, _j, _m, _2); \
	Two_Product_2Presplit(a1, a1hi, a1lo, b1, bhi, blo, _j, _0); \
	Two_Sum(_i, _0, _n, _0); \
	Two_Sum(_1, _0, _i, x2); \
	Two_Sum(_2, _i, _k, _1); \
	Two_Sum(_m, _k, _l, _2); \
	Two_Sum(_j, _n, _k, _0); \
	Two_Sum(_1, _0, _j, x3); \
	Two_Sum(_2, _j, _i, _1); \
	Two_Sum(_l, _i, _m, _2); \
	Two_Sum(_1, _k, _i, x4); \
	Two_Sum(_2, _i, _k, x5); \
	Two_Sum(_m, _k, x7, x6)
	*/

	/*
	#define Two_Square(a1, a0, x5, x4, x3, x2, x1, x0) \
	Square(a0, _j, x0); \
	_0 = a0 + a0; \
	Two_Product(a1, _0, _k, _1); \
	Two_One_Sum(_k, _1, _j, _l, _2, x1); \
	Square(a1, _j, _1); \
	Two_Two_Sum(_j, _1, _l, _2, x5, x4, x3, x2)
	*/

		//double splitter;     /* = 2^ceiling(p / 2) + 1.  Used to split floats in half. */
		//double epsilon;                /* = 2^(-p).  Used to estimate roundoff errors. */
									/* A set of coefficients used to calculate maximum roundoff errors.          */
		double resulterrbound = (3.0 + 8.0 * epsilon) * epsilon;
		double o3derrboundA = (7.0 + 56.0 * epsilon) * epsilon;
		double o3derrboundB = (3.0 + 28.0 * epsilon) * epsilon;
		double o3derrboundC = (26.0 + 288.0 * epsilon) * epsilon * epsilon;
		double isperrboundA = (16.0 + 224.0 * epsilon) * epsilon;
		double isperrboundB = (5.0 + 72.0 * epsilon) * epsilon;
		double isperrboundC = (71.0 + 1408.0 * epsilon) * epsilon * epsilon;


		/*double doublerand()
		{
			double result;
			double expo;
			long a, b, c;
			long i;

			a = random();
			b = random();
			c = random();
			result = (double)(a - 1073741824) * 8388608.0 + (double)(b >> 8);
			for (i = 512, expo = 2; i <= 131072; i *= 2, expo = expo * expo) {
				if (c & i) {
					result *= expo;
				}
			}
			return result;
		}

		double narrowdoublerand()
		{
			double result;
			double expo;
			long a, b, c;
			long i;

			a = random();
			b = random();
			c = random();
			result = (double)(a - 1073741824) * 8388608.0 + (double)(b >> 8);
			for (i = 512, expo = 2; i <= 2048; i *= 2, expo = expo * expo) {
				if (c & i) {
					result *= expo;
				}
			}
			return result;
		}

		double uniformdoublerand()
		{
			double result;
			long a, b;

			a = random();
			b = random();
			result = (double)(a - 1073741824) * 8388608.0 + (double)(b >> 8);
			return result;
		}

		float floatrand()
		{
			float result;
			float expo;
			long a, c;
			long i;

			a = random();
			c = random();
			result = (float)((a - 1073741824) >> 6);
			for (i = 512, expo = 2; i <= 16384; i *= 2, expo = expo * expo) {
				if (c & i) {
					result *= expo;
				}
			}
			return result;
		}

		float narrowfloatrand()
		{
			float result;
			float expo;
			long a, c;
			long i;

			a = random();
			c = random();
			result = (float)((a - 1073741824) >> 6);
			for (i = 512, expo = 2; i <= 2048; i *= 2, expo = expo * expo) {
				if (c & i) {
					result *= expo;
				}
			}
			return result;
		}

		float uniformfloatrand()
		{
			float result;
			long a;

			a = random();
			result = (float)((a - 1073741824) >> 6);
			return result;
		}
		*/

		int fast_expansion_sum_zeroelim(int elen, double *e, int flen, double *f, double *h)  /* h cannot be e or f. */
		{
			double Q;
			double Qnew;
			double hh;
			double bvirt;
			double avirt, bround, around;
			int eindex, findex, hindex;
			double enow, fnow;

			enow = e[0];
			fnow = f[0];
			eindex = findex = 0;
			if ((fnow > enow) == (fnow > -enow)) {
				Q = enow;
				enow = e[++eindex];
			}
			else {
				Q = fnow;
				fnow = f[++findex];
			}
			hindex = 0;
			if ((eindex < elen) && (findex < flen)) {
				if ((fnow > enow) == (fnow > -enow)) {
					Fast_Two_Sum(enow, Q, Qnew, hh);
					enow = e[++eindex];
				}
				else {
					Fast_Two_Sum(fnow, Q, Qnew, hh);
					fnow = f[++findex];
				}
				Q = Qnew;
				if (hh != 0.0) {
					h[hindex++] = hh;
				}
				while ((eindex < elen) && (findex < flen)) {
					if ((fnow > enow) == (fnow > -enow)) {
						Two_Sum(Q, enow, Qnew, hh);
						enow = e[++eindex];
					}
					else {
						Two_Sum(Q, fnow, Qnew, hh);
						fnow = f[++findex];
					}
					Q = Qnew;
					if (hh != 0.0) {
						h[hindex++] = hh;
					}
				}
			}
			while (eindex < elen) {
				Two_Sum(Q, enow, Qnew, hh);
				enow = e[++eindex];
				Q = Qnew;
				if (hh != 0.0) {
					h[hindex++] = hh;
				}
			}
			while (findex < flen) {
				Two_Sum(Q, fnow, Qnew, hh);
				fnow = f[++findex];
				Q = Qnew;
				if (hh != 0.0) {
					h[hindex++] = hh;
				}
			}
			if ((Q != 0.0) || (hindex == 0)) {
				h[hindex++] = Q;
			}
			return hindex;
		}


		int scale_expansion_zeroelim(int elen, double *e, double b, double *h)   /* e and h cannot be the same. */
		{
			double Q;
			double hh;
			//		double product1;
			//		double product0;
			int eindex, hindex;
			//		double enow;
			double bvirt;
			double avirt, bround, around;
			double c;
			double abig;
			double ahi, alo, bhi, blo;
			double err1, err2, err3;

			Split(b, bhi, blo);
			Two_Product_Presplit(e[0], b, bhi, blo, Q, hh);
			hindex = 0;
			if (hh != 0) {
				h[hindex++] = hh;
			}
			for (eindex = 1; eindex < elen; eindex++) {
			const double enow = e[eindex];
				double product1;
				double product0;
				Two_Product_Presplit(enow, b, bhi, blo, product1, product0);
				double sum;
				Two_Sum(Q, product0, sum, hh);
				if (hh != 0) {
					h[hindex++] = hh;
				}
				Fast_Two_Sum(product1, sum, Q, hh);
				if (hh != 0) {
					h[hindex++] = hh;
				}
			}
			if ((Q != 0.0) || (hindex == 0)) {
				h[hindex++] = Q;
			}
			return hindex;
		}

		double estimate(int elen, double *e)
		{
			double Q;
			int eindex;

			Q = e[0];
			for (eindex = 1; eindex < elen; eindex++) {
				Q += e[eindex];
			}
			return Q;
		}

	ATTRIBUTE_NO_SANITIZE_ADDRESS
		double orient3dadapt(std::array<Point3D, 4> const& points, double permanent)
		{
			double pa[3], pb[3], pc[3], pd[3];
			pa[0] = points[0].x;
			pa[1] = points[0].y;
			pa[2] = points[0].z;
			pb[0] = points[1].x;
			pb[1] = points[1].y;
			pb[2] = points[1].z;
			pc[0] = points[2].x;
			pc[1] = points[2].y;
			pc[2] = points[2].z;
			pd[0] = points[3].x;
			pd[1] = points[3].y;
			pd[2] = points[3].z;
			double adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz;
			double det, errbound;

			double bdxcdy1, cdxbdy1, cdxady1, adxcdy1, adxbdy1, bdxady1;
			double bdxcdy0, cdxbdy0, cdxady0, adxcdy0, adxbdy0, bdxady0;
			double bc[4], ca[4], ab[4];
			double bc3, ca3, ab3;
			double adet[8], bdet[8], cdet[8];
			int alen, blen, clen;
			double abdet[16];
			int ablen;
			double *finnow, *finother, *finswap;
			double fin1[192], fin2[192];
			int finlength;

			double adxtail, bdxtail, cdxtail;
			double adytail, bdytail, cdytail;
			double adztail, bdztail, cdztail;
			double at_blarge, at_clarge;
			double bt_clarge, bt_alarge;
			double ct_alarge, ct_blarge;
			double at_b[4], at_c[4], bt_c[4], bt_a[4], ct_a[4], ct_b[4];
			int at_blen, at_clen, bt_clen, bt_alen, ct_alen, ct_blen;
			double bdxt_cdy1, cdxt_bdy1, cdxt_ady1;
			double adxt_cdy1, adxt_bdy1, bdxt_ady1;
			double bdxt_cdy0, cdxt_bdy0, cdxt_ady0;
			double adxt_cdy0, adxt_bdy0, bdxt_ady0;
			double bdyt_cdx1, cdyt_bdx1, cdyt_adx1;
			double adyt_cdx1, adyt_bdx1, bdyt_adx1;
			double bdyt_cdx0, cdyt_bdx0, cdyt_adx0;
			double adyt_cdx0, adyt_bdx0, bdyt_adx0;
			double bct[8], cat[8], abt[8];
			int bctlen, catlen, abtlen;
			double bdxt_cdyt1, cdxt_bdyt1, cdxt_adyt1;
			double adxt_cdyt1, adxt_bdyt1, bdxt_adyt1;
			double bdxt_cdyt0, cdxt_bdyt0, cdxt_adyt0;
			double adxt_cdyt0, adxt_bdyt0, bdxt_adyt0;
			double u[4], v[12], w[16];
			double u3;
			int vlength, wlength;
			double negate;

			double bvirt;
			double avirt, bround, around;
			double c;
			double abig;
			double ahi, alo, bhi, blo;
			double err1, err2, err3;
			double _i, _j, _k;
			double _0;

			adx = static_cast<double>(pa[0] - pd[0]);
			bdx = static_cast<double>(pb[0] - pd[0]);
			cdx = static_cast<double>(pc[0] - pd[0]);
			ady = static_cast<double>(pa[1] - pd[1]);
			bdy = static_cast<double>(pb[1] - pd[1]);
			cdy = static_cast<double>(pc[1] - pd[1]);
			adz = static_cast<double>(pa[2] - pd[2]);
			bdz = static_cast<double>(pb[2] - pd[2]);
			cdz = static_cast<double>(pc[2] - pd[2]);

			Two_Product(bdx, cdy, bdxcdy1, bdxcdy0);
			Two_Product(cdx, bdy, cdxbdy1, cdxbdy0);
			Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc3, bc[2], bc[1], bc[0]);
			bc[3] = bc3;
			alen = scale_expansion_zeroelim(4, bc, adz, adet);

			Two_Product(cdx, ady, cdxady1, cdxady0);
			Two_Product(adx, cdy, adxcdy1, adxcdy0);
			Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca3, ca[2], ca[1], ca[0]);
			ca[3] = ca3;
			blen = scale_expansion_zeroelim(4, ca, bdz, bdet);

			Two_Product(adx, bdy, adxbdy1, adxbdy0);
			Two_Product(bdx, ady, bdxady1, bdxady0);
			Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab3, ab[2], ab[1], ab[0]);
			ab[3] = ab3;
			clen = scale_expansion_zeroelim(4, ab, cdz, cdet);

			ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
			finlength = fast_expansion_sum_zeroelim(ablen, abdet, clen, cdet, fin1);

			det = estimate(finlength, fin1);
			errbound = o3derrboundB * permanent;
			if ((det >= errbound) || (-det >= errbound)) {
				return det;
			}

			Two_Diff_Tail(pa[0], pd[0], adx, adxtail);
			Two_Diff_Tail(pb[0], pd[0], bdx, bdxtail);
			Two_Diff_Tail(pc[0], pd[0], cdx, cdxtail);
			Two_Diff_Tail(pa[1], pd[1], ady, adytail);
			Two_Diff_Tail(pb[1], pd[1], bdy, bdytail);
			Two_Diff_Tail(pc[1], pd[1], cdy, cdytail);
			Two_Diff_Tail(pa[2], pd[2], adz, adztail);
			Two_Diff_Tail(pb[2], pd[2], bdz, bdztail);
			Two_Diff_Tail(pc[2], pd[2], cdz, cdztail);

			if ((adxtail == 0.0) && (bdxtail == 0.0) && (cdxtail == 0.0)
				&& (adytail == 0.0) && (bdytail == 0.0) && (cdytail == 0.0)
				&& (adztail == 0.0) && (bdztail == 0.0) && (cdztail == 0.0)) {
				return det;
			}

			errbound = o3derrboundC * permanent + resulterrbound * Absolute(det);
			det += (adz * ((bdx * cdytail + cdy * bdxtail)
				- (bdy * cdxtail + cdx * bdytail))
				+ adztail * (bdx * cdy - bdy * cdx))
				+ (bdz * ((cdx * adytail + ady * cdxtail)
					- (cdy * adxtail + adx * cdytail))
					+ bdztail * (cdx * ady - cdy * adx))
				+ (cdz * ((adx * bdytail + bdy * adxtail)
					- (ady * bdxtail + bdx * adytail))
					+ cdztail * (adx * bdy - ady * bdx));
			if ((det >= errbound) || (-det >= errbound)) {
				return det;
			}

			finnow = fin1;
			finother = fin2;

			if (adxtail == 0.0) {
				if (adytail == 0.0) {
					at_b[0] = 0.0;
					at_blen = 1;
					at_c[0] = 0.0;
					at_clen = 1;
				}
				else {
					negate = -adytail;
					Two_Product(negate, bdx, at_blarge, at_b[0]);
					at_b[1] = at_blarge;
					at_blen = 2;
					Two_Product(adytail, cdx, at_clarge, at_c[0]);
					at_c[1] = at_clarge;
					at_clen = 2;
				}
			}
			else {
				if (adytail == 0.0) {
					Two_Product(adxtail, bdy, at_blarge, at_b[0]);
					at_b[1] = at_blarge;
					at_blen = 2;
					negate = -adxtail;
					Two_Product(negate, cdy, at_clarge, at_c[0]);
					at_c[1] = at_clarge;
					at_clen = 2;
				}
				else {
					Two_Product(adxtail, bdy, adxt_bdy1, adxt_bdy0);
					Two_Product(adytail, bdx, adyt_bdx1, adyt_bdx0);
					Two_Two_Diff(adxt_bdy1, adxt_bdy0, adyt_bdx1, adyt_bdx0,
						at_blarge, at_b[2], at_b[1], at_b[0]);
					at_b[3] = at_blarge;
					at_blen = 4;
					Two_Product(adytail, cdx, adyt_cdx1, adyt_cdx0);
					Two_Product(adxtail, cdy, adxt_cdy1, adxt_cdy0);
					Two_Two_Diff(adyt_cdx1, adyt_cdx0, adxt_cdy1, adxt_cdy0,
						at_clarge, at_c[2], at_c[1], at_c[0]);
					at_c[3] = at_clarge;
					at_clen = 4;
				}
			}
			if (bdxtail == 0.0) {
				if (bdytail == 0.0) {
					bt_c[0] = 0.0;
					bt_clen = 1;
					bt_a[0] = 0.0;
					bt_alen = 1;
				}
				else {
					negate = -bdytail;
					Two_Product(negate, cdx, bt_clarge, bt_c[0]);
					bt_c[1] = bt_clarge;
					bt_clen = 2;
					Two_Product(bdytail, adx, bt_alarge, bt_a[0]);
					bt_a[1] = bt_alarge;
					bt_alen = 2;
				}
			}
			else {
				if (bdytail == 0.0) {
					Two_Product(bdxtail, cdy, bt_clarge, bt_c[0]);
					bt_c[1] = bt_clarge;
					bt_clen = 2;
					negate = -bdxtail;
					Two_Product(negate, ady, bt_alarge, bt_a[0]);
					bt_a[1] = bt_alarge;
					bt_alen = 2;
				}
				else {
					Two_Product(bdxtail, cdy, bdxt_cdy1, bdxt_cdy0);
					Two_Product(bdytail, cdx, bdyt_cdx1, bdyt_cdx0);
					Two_Two_Diff(bdxt_cdy1, bdxt_cdy0, bdyt_cdx1, bdyt_cdx0,
						bt_clarge, bt_c[2], bt_c[1], bt_c[0]);
					bt_c[3] = bt_clarge;
					bt_clen = 4;
					Two_Product(bdytail, adx, bdyt_adx1, bdyt_adx0);
					Two_Product(bdxtail, ady, bdxt_ady1, bdxt_ady0);
					Two_Two_Diff(bdyt_adx1, bdyt_adx0, bdxt_ady1, bdxt_ady0,
						bt_alarge, bt_a[2], bt_a[1], bt_a[0]);
					bt_a[3] = bt_alarge;
					bt_alen = 4;
				}
			}
			if (cdxtail == 0.0) {
				if (cdytail == 0.0) {
					ct_a[0] = 0.0;
					ct_alen = 1;
					ct_b[0] = 0.0;
					ct_blen = 1;
				}
				else {
					negate = -cdytail;
					Two_Product(negate, adx, ct_alarge, ct_a[0]);
					ct_a[1] = ct_alarge;
					ct_alen = 2;
					Two_Product(cdytail, bdx, ct_blarge, ct_b[0]);
					ct_b[1] = ct_blarge;
					ct_blen = 2;
				}
			}
			else {
				if (cdytail == 0.0) {
					Two_Product(cdxtail, ady, ct_alarge, ct_a[0]);
					ct_a[1] = ct_alarge;
					ct_alen = 2;
					negate = -cdxtail;
					Two_Product(negate, bdy, ct_blarge, ct_b[0]);
					ct_b[1] = ct_blarge;
					ct_blen = 2;
				}
				else {
					Two_Product(cdxtail, ady, cdxt_ady1, cdxt_ady0);
					Two_Product(cdytail, adx, cdyt_adx1, cdyt_adx0);
					Two_Two_Diff(cdxt_ady1, cdxt_ady0, cdyt_adx1, cdyt_adx0,
						ct_alarge, ct_a[2], ct_a[1], ct_a[0]);
					ct_a[3] = ct_alarge;
					ct_alen = 4;
					Two_Product(cdytail, bdx, cdyt_bdx1, cdyt_bdx0);
					Two_Product(cdxtail, bdy, cdxt_bdy1, cdxt_bdy0);
					Two_Two_Diff(cdyt_bdx1, cdyt_bdx0, cdxt_bdy1, cdxt_bdy0,
						ct_blarge, ct_b[2], ct_b[1], ct_b[0]);
					ct_b[3] = ct_blarge;
					ct_blen = 4;
				}
			}

			bctlen = fast_expansion_sum_zeroelim(bt_clen, bt_c, ct_blen, ct_b, bct);
			wlength = scale_expansion_zeroelim(bctlen, bct, adz, w);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
				finother);
			finswap = finnow; finnow = finother; finother = finswap;

			catlen = fast_expansion_sum_zeroelim(ct_alen, ct_a, at_clen, at_c, cat);
			wlength = scale_expansion_zeroelim(catlen, cat, bdz, w);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
				finother);
			finswap = finnow; finnow = finother; finother = finswap;

			abtlen = fast_expansion_sum_zeroelim(at_blen, at_b, bt_alen, bt_a, abt);
			wlength = scale_expansion_zeroelim(abtlen, abt, cdz, w);
			finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
				finother);
			finswap = finnow; finnow = finother; finother = finswap;

			if (adztail != 0.0) {
				vlength = scale_expansion_zeroelim(4, bc, adztail, v);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
					finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
			if (bdztail != 0.0) {
				vlength = scale_expansion_zeroelim(4, ca, bdztail, v);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
					finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
			if (cdztail != 0.0) {
				vlength = scale_expansion_zeroelim(4, ab, cdztail, v);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, vlength, v,
					finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}

			if (adxtail != 0.0) {
				if (bdytail != 0.0) {
					Two_Product(adxtail, bdytail, adxt_bdyt1, adxt_bdyt0);
					Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdz, u3, u[2], u[1], u[0]);
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
						finother);
					finswap = finnow; finnow = finother; finother = finswap;
					if (cdztail != 0.0) {
						Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdztail, u3, u[2], u[1], u[0]);
						u[3] = u3;
						finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
							finother);
						finswap = finnow; finnow = finother; finother = finswap;
					}
				}
				if (cdytail != 0.0) {
					negate = -adxtail;
					Two_Product(negate, cdytail, adxt_cdyt1, adxt_cdyt0);
					Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdz, u3, u[2], u[1], u[0]);
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
						finother);
					finswap = finnow; finnow = finother; finother = finswap;
					if (bdztail != 0.0) {
						Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdztail, u3, u[2], u[1], u[0]);
						u[3] = u3;
						finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
							finother);
						finswap = finnow; finnow = finother; finother = finswap;
					}
				}
			}
			if (bdxtail != 0.0) {
				if (cdytail != 0.0) {
					Two_Product(bdxtail, cdytail, bdxt_cdyt1, bdxt_cdyt0);
					Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adz, u3, u[2], u[1], u[0]);
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
						finother);
					finswap = finnow; finnow = finother; finother = finswap;
					if (adztail != 0.0) {
						Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adztail, u3, u[2], u[1], u[0]);
						u[3] = u3;
						finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
							finother);
						finswap = finnow; finnow = finother; finother = finswap;
					}
				}
				if (adytail != 0.0) {
					negate = -bdxtail;
					Two_Product(negate, adytail, bdxt_adyt1, bdxt_adyt0);
					Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdz, u3, u[2], u[1], u[0]);
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
						finother);
					finswap = finnow; finnow = finother; finother = finswap;
					if (cdztail != 0.0) {
						Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdztail, u3, u[2], u[1], u[0]);
						u[3] = u3;
						finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
							finother);
						finswap = finnow; finnow = finother; finother = finswap;
					}
				}
			}
			if (cdxtail != 0.0) {
				if (adytail != 0.0) {
					Two_Product(cdxtail, adytail, cdxt_adyt1, cdxt_adyt0);
					Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdz, u3, u[2], u[1], u[0]);
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
						finother);
					finswap = finnow; finnow = finother; finother = finswap;
					if (bdztail != 0.0) {
						Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdztail, u3, u[2], u[1], u[0]);
						u[3] = u3;
						finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
							finother);
						finswap = finnow; finnow = finother; finother = finswap;
					}
				}
				if (bdytail != 0.0) {
					negate = -cdxtail;
					Two_Product(negate, bdytail, cdxt_bdyt1, cdxt_bdyt0);
					Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adz, u3, u[2], u[1], u[0]);
					u[3] = u3;
					finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
						finother);
					finswap = finnow; finnow = finother; finother = finswap;
					if (adztail != 0.0) {
						Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adztail, u3, u[2], u[1], u[0]);
						u[3] = u3;
						finlength = fast_expansion_sum_zeroelim(finlength, finnow, 4, u,
							finother);
						finswap = finnow; finnow = finother; finother = finswap;
					}
				}
			}

			if (adztail != 0.0) {
				wlength = scale_expansion_zeroelim(bctlen, bct, adztail, w);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
					finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
			if (bdztail != 0.0) {
				wlength = scale_expansion_zeroelim(catlen, cat, bdztail, w);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
					finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}
			if (cdztail != 0.0) {
				wlength = scale_expansion_zeroelim(abtlen, abt, cdztail, w);
				finlength = fast_expansion_sum_zeroelim(finlength, finnow, wlength, w,
					finother);
				finswap = finnow; finnow = finother; finother = finswap;
			}

			return finnow[finlength - 1];
		}

		

		double insphereexact(double *pa, double * pb, double * pc, double * pd, double * pe)
		{
			double axby1, bxcy1, cxdy1, dxey1, exay1;
			double bxay1, cxby1, dxcy1, exdy1, axey1;
			double axcy1, bxdy1, cxey1, dxay1, exby1;
			double cxay1, dxby1, excy1, axdy1, bxey1;
			double axby0, bxcy0, cxdy0, dxey0, exay0;
			double bxay0, cxby0, dxcy0, exdy0, axey0;
			double axcy0, bxdy0, cxey0, dxay0, exby0;
			double cxay0, dxby0, excy0, axdy0, bxey0;
			double ab[4], bc[4], cd[4], de[4], ea[4];
			double ac[4], bd[4], ce[4], da[4], eb[4];
			double temp8a[8], temp8b[8], temp16[16];
			int temp8alen, temp8blen, temp16len;
			double abc[24], bcd[24], cde[24], dea[24], eab[24];
			double abd[24], bce[24], cda[24], deb[24], eac[24];
			int abclen, bcdlen, cdelen, dealen, eablen;
			int abdlen, bcelen, cdalen, deblen, eaclen;
			double temp48a[48], temp48b[48];
			int temp48alen, temp48blen;
			double abcd[96], bcde[96], cdea[96], deab[96], eabc[96];
			int abcdlen, bcdelen, cdealen, deablen, eabclen;
			double temp192[192];
			double det384x[384], det384y[384], det384z[384];
			int xlen, ylen, zlen;
			double detxy[768];
			int xylen;
			double adet[1152], bdet[1152], cdet[1152], ddet[1152], edet[1152];
			int alen, blen, clen, dlen, elen;
			double abdet[2304], cddet[2304], cdedet[3456];
			int ablen, cdlen;
			double deter[5760];
			int deterlen;
			int i;

			double bvirt;
			double avirt, bround, around;
			double c;
			double abig;
			double ahi, alo, bhi, blo;
			double err1, err2, err3;
			double _i, _j;
			double _0;

			Two_Product(pa[0], pb[1], axby1, axby0);
			Two_Product(pb[0], pa[1], bxay1, bxay0);
			Two_Two_Diff(axby1, axby0, bxay1, bxay0, ab[3], ab[2], ab[1], ab[0]);

			Two_Product(pb[0], pc[1], bxcy1, bxcy0);
			Two_Product(pc[0], pb[1], cxby1, cxby0);
			Two_Two_Diff(bxcy1, bxcy0, cxby1, cxby0, bc[3], bc[2], bc[1], bc[0]);

			Two_Product(pc[0], pd[1], cxdy1, cxdy0);
			Two_Product(pd[0], pc[1], dxcy1, dxcy0);
			Two_Two_Diff(cxdy1, cxdy0, dxcy1, dxcy0, cd[3], cd[2], cd[1], cd[0]);

			Two_Product(pd[0], pe[1], dxey1, dxey0);
			Two_Product(pe[0], pd[1], exdy1, exdy0);
			Two_Two_Diff(dxey1, dxey0, exdy1, exdy0, de[3], de[2], de[1], de[0]);

			Two_Product(pe[0], pa[1], exay1, exay0);
			Two_Product(pa[0], pe[1], axey1, axey0);
			Two_Two_Diff(exay1, exay0, axey1, axey0, ea[3], ea[2], ea[1], ea[0]);

			Two_Product(pa[0], pc[1], axcy1, axcy0);
			Two_Product(pc[0], pa[1], cxay1, cxay0);
			Two_Two_Diff(axcy1, axcy0, cxay1, cxay0, ac[3], ac[2], ac[1], ac[0]);

			Two_Product(pb[0], pd[1], bxdy1, bxdy0);
			Two_Product(pd[0], pb[1], dxby1, dxby0);
			Two_Two_Diff(bxdy1, bxdy0, dxby1, dxby0, bd[3], bd[2], bd[1], bd[0]);

			Two_Product(pc[0], pe[1], cxey1, cxey0);
			Two_Product(pe[0], pc[1], excy1, excy0);
			Two_Two_Diff(cxey1, cxey0, excy1, excy0, ce[3], ce[2], ce[1], ce[0]);

			Two_Product(pd[0], pa[1], dxay1, dxay0);
			Two_Product(pa[0], pd[1], axdy1, axdy0);
			Two_Two_Diff(dxay1, dxay0, axdy1, axdy0, da[3], da[2], da[1], da[0]);

			Two_Product(pe[0], pb[1], exby1, exby0);
			Two_Product(pb[0], pe[1], bxey1, bxey0);
			Two_Two_Diff(exby1, exby0, bxey1, bxey0, eb[3], eb[2], eb[1], eb[0]);

			temp8alen = scale_expansion_zeroelim(4, bc, pa[2], temp8a);
			temp8blen = scale_expansion_zeroelim(4, ac, -pb[2], temp8b);
			temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
			temp8alen = scale_expansion_zeroelim(4, ab, pc[2], temp8a);
			abclen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				abc);

			temp8alen = scale_expansion_zeroelim(4, cd, pb[2], temp8a);
			temp8blen = scale_expansion_zeroelim(4, bd, -pc[2], temp8b);
			temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
			temp8alen = scale_expansion_zeroelim(4, bc, pd[2], temp8a);
			bcdlen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				bcd);

			temp8alen = scale_expansion_zeroelim(4, de, pc[2], temp8a);
			temp8blen = scale_expansion_zeroelim(4, ce, -pd[2], temp8b);
			temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
			temp8alen = scale_expansion_zeroelim(4, cd, pe[2], temp8a);
			cdelen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				cde);

			temp8alen = scale_expansion_zeroelim(4, ea, pd[2], temp8a);
			temp8blen = scale_expansion_zeroelim(4, da, -pe[2], temp8b);
			temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
			temp8alen = scale_expansion_zeroelim(4, de, pa[2], temp8a);
			dealen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				dea);

			temp8alen = scale_expansion_zeroelim(4, ab, pe[2], temp8a);
			temp8blen = scale_expansion_zeroelim(4, eb, -pa[2], temp8b);
			temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
			temp8alen = scale_expansion_zeroelim(4, ea, pb[2], temp8a);
			eablen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				eab);

			temp8alen = scale_expansion_zeroelim(4, bd, pa[2], temp8a);
			temp8blen = scale_expansion_zeroelim(4, da, pb[2], temp8b);
			temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
			temp8alen = scale_expansion_zeroelim(4, ab, pd[2], temp8a);
			abdlen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				abd);

			temp8alen = scale_expansion_zeroelim(4, ce, pb[2], temp8a);
			temp8blen = scale_expansion_zeroelim(4, eb, pc[2], temp8b);
			temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
			temp8alen = scale_expansion_zeroelim(4, bc, pe[2], temp8a);
			bcelen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				bce);

			temp8alen = scale_expansion_zeroelim(4, da, pc[2], temp8a);
			temp8blen = scale_expansion_zeroelim(4, ac, pd[2], temp8b);
			temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
			temp8alen = scale_expansion_zeroelim(4, cd, pa[2], temp8a);
			cdalen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				cda);

			temp8alen = scale_expansion_zeroelim(4, eb, pd[2], temp8a);
			temp8blen = scale_expansion_zeroelim(4, bd, pe[2], temp8b);
			temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
			temp8alen = scale_expansion_zeroelim(4, de, pb[2], temp8a);
			deblen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				deb);

			temp8alen = scale_expansion_zeroelim(4, ac, pe[2], temp8a);
			temp8blen = scale_expansion_zeroelim(4, ce, pa[2], temp8b);
			temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp8blen, temp8b,
				temp16);
			temp8alen = scale_expansion_zeroelim(4, ea, pc[2], temp8a);
			eaclen = fast_expansion_sum_zeroelim(temp8alen, temp8a, temp16len, temp16,
				eac);

			temp48alen = fast_expansion_sum_zeroelim(cdelen, cde, bcelen, bce, temp48a);
			temp48blen = fast_expansion_sum_zeroelim(deblen, deb, bcdlen, bcd, temp48b);
			for (i = 0; i < temp48blen; i++) {
				temp48b[i] = -temp48b[i];
			}
			bcdelen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
				temp48blen, temp48b, bcde);
			xlen = scale_expansion_zeroelim(bcdelen, bcde, pa[0], temp192);
			xlen = scale_expansion_zeroelim(xlen, temp192, pa[0], det384x);
			ylen = scale_expansion_zeroelim(bcdelen, bcde, pa[1], temp192);
			ylen = scale_expansion_zeroelim(ylen, temp192, pa[1], det384y);
			zlen = scale_expansion_zeroelim(bcdelen, bcde, pa[2], temp192);
			zlen = scale_expansion_zeroelim(zlen, temp192, pa[2], det384z);
			xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
			alen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, adet);

			temp48alen = fast_expansion_sum_zeroelim(dealen, dea, cdalen, cda, temp48a);
			temp48blen = fast_expansion_sum_zeroelim(eaclen, eac, cdelen, cde, temp48b);
			for (i = 0; i < temp48blen; i++) {
				temp48b[i] = -temp48b[i];
			}
			cdealen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
				temp48blen, temp48b, cdea);
			xlen = scale_expansion_zeroelim(cdealen, cdea, pb[0], temp192);
			xlen = scale_expansion_zeroelim(xlen, temp192, pb[0], det384x);
			ylen = scale_expansion_zeroelim(cdealen, cdea, pb[1], temp192);
			ylen = scale_expansion_zeroelim(ylen, temp192, pb[1], det384y);
			zlen = scale_expansion_zeroelim(cdealen, cdea, pb[2], temp192);
			zlen = scale_expansion_zeroelim(zlen, temp192, pb[2], det384z);
			xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
			blen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, bdet);

			temp48alen = fast_expansion_sum_zeroelim(eablen, eab, deblen, deb, temp48a);
			temp48blen = fast_expansion_sum_zeroelim(abdlen, abd, dealen, dea, temp48b);
			for (i = 0; i < temp48blen; i++) {
				temp48b[i] = -temp48b[i];
			}
			deablen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
				temp48blen, temp48b, deab);
			xlen = scale_expansion_zeroelim(deablen, deab, pc[0], temp192);
			xlen = scale_expansion_zeroelim(xlen, temp192, pc[0], det384x);
			ylen = scale_expansion_zeroelim(deablen, deab, pc[1], temp192);
			ylen = scale_expansion_zeroelim(ylen, temp192, pc[1], det384y);
			zlen = scale_expansion_zeroelim(deablen, deab, pc[2], temp192);
			zlen = scale_expansion_zeroelim(zlen, temp192, pc[2], det384z);
			xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
			clen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, cdet);

			temp48alen = fast_expansion_sum_zeroelim(abclen, abc, eaclen, eac, temp48a);
			temp48blen = fast_expansion_sum_zeroelim(bcelen, bce, eablen, eab, temp48b);
			for (i = 0; i < temp48blen; i++) {
				temp48b[i] = -temp48b[i];
			}
			eabclen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
				temp48blen, temp48b, eabc);
			xlen = scale_expansion_zeroelim(eabclen, eabc, pd[0], temp192);
			xlen = scale_expansion_zeroelim(xlen, temp192, pd[0], det384x);
			ylen = scale_expansion_zeroelim(eabclen, eabc, pd[1], temp192);
			ylen = scale_expansion_zeroelim(ylen, temp192, pd[1], det384y);
			zlen = scale_expansion_zeroelim(eabclen, eabc, pd[2], temp192);
			zlen = scale_expansion_zeroelim(zlen, temp192, pd[2], det384z);
			xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
			dlen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, ddet);

			temp48alen = fast_expansion_sum_zeroelim(bcdlen, bcd, abdlen, abd, temp48a);
			temp48blen = fast_expansion_sum_zeroelim(cdalen, cda, abclen, abc, temp48b);
			for (i = 0; i < temp48blen; i++) {
				temp48b[i] = -temp48b[i];
			}
			abcdlen = fast_expansion_sum_zeroelim(temp48alen, temp48a,
				temp48blen, temp48b, abcd);
			xlen = scale_expansion_zeroelim(abcdlen, abcd, pe[0], temp192);
			xlen = scale_expansion_zeroelim(xlen, temp192, pe[0], det384x);
			ylen = scale_expansion_zeroelim(abcdlen, abcd, pe[1], temp192);
			ylen = scale_expansion_zeroelim(ylen, temp192, pe[1], det384y);
			zlen = scale_expansion_zeroelim(abcdlen, abcd, pe[2], temp192);
			zlen = scale_expansion_zeroelim(zlen, temp192, pe[2], det384z);
			xylen = fast_expansion_sum_zeroelim(xlen, det384x, ylen, det384y, detxy);
			elen = fast_expansion_sum_zeroelim(xylen, detxy, zlen, det384z, edet);

			ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
			cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
			cdelen = fast_expansion_sum_zeroelim(cdlen, cddet, elen, edet, cdedet);
			deterlen = fast_expansion_sum_zeroelim(ablen, abdet, cdelen, cdedet, deter);

			return deter[deterlen - 1];
		}


		double insphereadapt(double *pa, double * pb, double * pc, double * pd, double * pe, double permanent)
		{
			double aex, bex, cex, dex, aey, bey, cey, dey, aez, bez, cez, dez;
			double det, errbound;

			double aexbey1, bexaey1, bexcey1, cexbey1;
			double cexdey1, dexcey1, dexaey1, aexdey1;
			double aexcey1, cexaey1, bexdey1, dexbey1;
			double aexbey0, bexaey0, bexcey0, cexbey0;
			double cexdey0, dexcey0, dexaey0, aexdey0;
			double aexcey0, cexaey0, bexdey0, dexbey0;
			double ab[4], bc[4], cd[4], da[4], ac[4], bd[4];
			double ab3, bc3, cd3, da3, ac3, bd3;
			double abeps, bceps, cdeps, daeps, aceps, bdeps;
			double temp8a[8], temp8b[8], temp8c[8], temp16[16], temp24[24], temp48[48];
			int temp8alen, temp8blen, temp8clen, temp16len, temp24len, temp48len;
			double xdet[96], ydet[96], zdet[96], xydet[192];
			int xlen, ylen, zlen, xylen;
			double adet[288], bdet[288], cdet[288], ddet[288];
			int alen, blen, clen, dlen;
			double abdet[576], cddet[576];
			int ablen, cdlen;
			double fin1[1152];
			int finlength;

			double aextail, bextail, cextail, dextail;
			double aeytail, beytail, ceytail, deytail;
			double aeztail, beztail, ceztail, deztail;

			double bvirt;
			double avirt, bround, around;
			double c;
			double abig;
			double ahi, alo, bhi, blo;
			double err1, err2, err3;
			double _i, _j;
			double _0;

			aex = static_cast<double>(pa[0] - pe[0]);
			bex = static_cast<double>(pb[0] - pe[0]);
			cex = static_cast<double>(pc[0] - pe[0]);
			dex = static_cast<double>(pd[0] - pe[0]);
			aey = static_cast<double>(pa[1] - pe[1]);
			bey = static_cast<double>(pb[1] - pe[1]);
			cey = static_cast<double>(pc[1] - pe[1]);
			dey = static_cast<double>(pd[1] - pe[1]);
			aez = static_cast<double>(pa[2] - pe[2]);
			bez = static_cast<double>(pb[2] - pe[2]);
			cez = static_cast<double>(pc[2] - pe[2]);
			dez = static_cast<double>(pd[2] - pe[2]);

			Two_Product(aex, bey, aexbey1, aexbey0);
			Two_Product(bex, aey, bexaey1, bexaey0);
			Two_Two_Diff(aexbey1, aexbey0, bexaey1, bexaey0, ab3, ab[2], ab[1], ab[0]);
			ab[3] = ab3;

			Two_Product(bex, cey, bexcey1, bexcey0);
			Two_Product(cex, bey, cexbey1, cexbey0);
			Two_Two_Diff(bexcey1, bexcey0, cexbey1, cexbey0, bc3, bc[2], bc[1], bc[0]);
			bc[3] = bc3;

			Two_Product(cex, dey, cexdey1, cexdey0);
			Two_Product(dex, cey, dexcey1, dexcey0);
			Two_Two_Diff(cexdey1, cexdey0, dexcey1, dexcey0, cd3, cd[2], cd[1], cd[0]);
			cd[3] = cd3;

			Two_Product(dex, aey, dexaey1, dexaey0);
			Two_Product(aex, dey, aexdey1, aexdey0);
			Two_Two_Diff(dexaey1, dexaey0, aexdey1, aexdey0, da3, da[2], da[1], da[0]);
			da[3] = da3;

			Two_Product(aex, cey, aexcey1, aexcey0);
			Two_Product(cex, aey, cexaey1, cexaey0);
			Two_Two_Diff(aexcey1, aexcey0, cexaey1, cexaey0, ac3, ac[2], ac[1], ac[0]);
			ac[3] = ac3;

			Two_Product(bex, dey, bexdey1, bexdey0);
			Two_Product(dex, bey, dexbey1, dexbey0);
			Two_Two_Diff(bexdey1, bexdey0, dexbey1, dexbey0, bd3, bd[2], bd[1], bd[0]);
			bd[3] = bd3;

			temp8alen = scale_expansion_zeroelim(4, cd, bez, temp8a);
			temp8blen = scale_expansion_zeroelim(4, bd, -cez, temp8b);
			temp8clen = scale_expansion_zeroelim(4, bc, dez, temp8c);
			temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
				temp8blen, temp8b, temp16);
			temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
				temp16len, temp16, temp24);
			temp48len = scale_expansion_zeroelim(temp24len, temp24, aex, temp48);
			xlen = scale_expansion_zeroelim(temp48len, temp48, -aex, xdet);
			temp48len = scale_expansion_zeroelim(temp24len, temp24, aey, temp48);
			ylen = scale_expansion_zeroelim(temp48len, temp48, -aey, ydet);
			temp48len = scale_expansion_zeroelim(temp24len, temp24, aez, temp48);
			zlen = scale_expansion_zeroelim(temp48len, temp48, -aez, zdet);
			xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
			alen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, adet);

			temp8alen = scale_expansion_zeroelim(4, da, cez, temp8a);
			temp8blen = scale_expansion_zeroelim(4, ac, dez, temp8b);
			temp8clen = scale_expansion_zeroelim(4, cd, aez, temp8c);
			temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
				temp8blen, temp8b, temp16);
			temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
				temp16len, temp16, temp24);
			temp48len = scale_expansion_zeroelim(temp24len, temp24, bex, temp48);
			xlen = scale_expansion_zeroelim(temp48len, temp48, bex, xdet);
			temp48len = scale_expansion_zeroelim(temp24len, temp24, bey, temp48);
			ylen = scale_expansion_zeroelim(temp48len, temp48, bey, ydet);
			temp48len = scale_expansion_zeroelim(temp24len, temp24, bez, temp48);
			zlen = scale_expansion_zeroelim(temp48len, temp48, bez, zdet);
			xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
			blen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, bdet);

			temp8alen = scale_expansion_zeroelim(4, ab, dez, temp8a);
			temp8blen = scale_expansion_zeroelim(4, bd, aez, temp8b);
			temp8clen = scale_expansion_zeroelim(4, da, bez, temp8c);
			temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
				temp8blen, temp8b, temp16);
			temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
				temp16len, temp16, temp24);
			temp48len = scale_expansion_zeroelim(temp24len, temp24, cex, temp48);
			xlen = scale_expansion_zeroelim(temp48len, temp48, -cex, xdet);
			temp48len = scale_expansion_zeroelim(temp24len, temp24, cey, temp48);
			ylen = scale_expansion_zeroelim(temp48len, temp48, -cey, ydet);
			temp48len = scale_expansion_zeroelim(temp24len, temp24, cez, temp48);
			zlen = scale_expansion_zeroelim(temp48len, temp48, -cez, zdet);
			xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
			clen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, cdet);

			temp8alen = scale_expansion_zeroelim(4, bc, aez, temp8a);
			temp8blen = scale_expansion_zeroelim(4, ac, -bez, temp8b);
			temp8clen = scale_expansion_zeroelim(4, ab, cez, temp8c);
			temp16len = fast_expansion_sum_zeroelim(temp8alen, temp8a,
				temp8blen, temp8b, temp16);
			temp24len = fast_expansion_sum_zeroelim(temp8clen, temp8c,
				temp16len, temp16, temp24);
			temp48len = scale_expansion_zeroelim(temp24len, temp24, dex, temp48);
			xlen = scale_expansion_zeroelim(temp48len, temp48, dex, xdet);
			temp48len = scale_expansion_zeroelim(temp24len, temp24, dey, temp48);
			ylen = scale_expansion_zeroelim(temp48len, temp48, dey, ydet);
			temp48len = scale_expansion_zeroelim(temp24len, temp24, dez, temp48);
			zlen = scale_expansion_zeroelim(temp48len, temp48, dez, zdet);
			xylen = fast_expansion_sum_zeroelim(xlen, xdet, ylen, ydet, xydet);
			dlen = fast_expansion_sum_zeroelim(xylen, xydet, zlen, zdet, ddet);

			ablen = fast_expansion_sum_zeroelim(alen, adet, blen, bdet, abdet);
			cdlen = fast_expansion_sum_zeroelim(clen, cdet, dlen, ddet, cddet);
			finlength = fast_expansion_sum_zeroelim(ablen, abdet, cdlen, cddet, fin1);

			det = estimate(finlength, fin1);
			errbound = isperrboundB * permanent;
			if ((det >= errbound) || (-det >= errbound)) {
				return det;
			}

			Two_Diff_Tail(pa[0], pe[0], aex, aextail);
			Two_Diff_Tail(pa[1], pe[1], aey, aeytail);
			Two_Diff_Tail(pa[2], pe[2], aez, aeztail);
			Two_Diff_Tail(pb[0], pe[0], bex, bextail);
			Two_Diff_Tail(pb[1], pe[1], bey, beytail);
			Two_Diff_Tail(pb[2], pe[2], bez, beztail);
			Two_Diff_Tail(pc[0], pe[0], cex, cextail);
			Two_Diff_Tail(pc[1], pe[1], cey, ceytail);
			Two_Diff_Tail(pc[2], pe[2], cez, ceztail);
			Two_Diff_Tail(pd[0], pe[0], dex, dextail);
			Two_Diff_Tail(pd[1], pe[1], dey, deytail);
			Two_Diff_Tail(pd[2], pe[2], dez, deztail);
			if ((aextail == 0.0) && (aeytail == 0.0) && (aeztail == 0.0)
				&& (bextail == 0.0) && (beytail == 0.0) && (beztail == 0.0)
				&& (cextail == 0.0) && (ceytail == 0.0) && (ceztail == 0.0)
				&& (dextail == 0.0) && (deytail == 0.0) && (deztail == 0.0)) {
				return det;
			}

			errbound = isperrboundC * permanent + resulterrbound * Absolute(det);
			abeps = (aex * beytail + bey * aextail)
				- (aey * bextail + bex * aeytail);
			bceps = (bex * ceytail + cey * bextail)
				- (bey * cextail + cex * beytail);
			cdeps = (cex * deytail + dey * cextail)
				- (cey * dextail + dex * ceytail);
			daeps = (dex * aeytail + aey * dextail)
				- (dey * aextail + aex * deytail);
			aceps = (aex * ceytail + cey * aextail)
				- (aey * cextail + cex * aeytail);
			bdeps = (bex * deytail + dey * bextail)
				- (bey * dextail + dex * beytail);
			det += (((bex * bex + bey * bey + bez * bez)
				* ((cez * daeps + dez * aceps + aez * cdeps)
					+ (ceztail * da3 + deztail * ac3 + aeztail * cd3))
				+ (dex * dex + dey * dey + dez * dez)
				* ((aez * bceps - bez * aceps + cez * abeps)
					+ (aeztail * bc3 - beztail * ac3 + ceztail * ab3)))
				- ((aex * aex + aey * aey + aez * aez)
					* ((bez * cdeps - cez * bdeps + dez * bceps)
						+ (beztail * cd3 - ceztail * bd3 + deztail * bc3))
					+ (cex * cex + cey * cey + cez * cez)
					* ((dez * abeps + aez * bdeps + bez * daeps)
						+ (deztail * ab3 + aeztail * bd3 + beztail * da3))))
				+ 2.0 * (((bex * bextail + bey * beytail + bez * beztail)
					* (cez * da3 + dez * ac3 + aez * cd3)
					+ (dex * dextail + dey * deytail + dez * deztail)
					* (aez * bc3 - bez * ac3 + cez * ab3))
					- ((aex * aextail + aey * aeytail + aez * aeztail)
						* (bez * cd3 - cez * bd3 + dez * bc3)
						+ (cex * cextail + cey * ceytail + cez * ceztail)
						* (dez * ab3 + aez * bd3 + bez * da3)));
			if ((det >= errbound) || (-det >= errbound)) {
				return det;
			}

			return insphereexact(pa, pb, pc, pd, pe);
		}

		double insphere(double *pa, double * pb, double * pc, double * pd, double * pe)
		{
			double aex, bex, cex, dex;
			double aey, bey, cey, dey;
			double aez, bez, cez, dez;
			double aexbey, bexaey, bexcey, cexbey, cexdey, dexcey, dexaey, aexdey;
			double aexcey, cexaey, bexdey, dexbey;
			double alift, blift, clift, dlift;
			double ab, bc, cd, da, ac, bd;
			double abc, bcd, cda, dab;
			double aezplus, bezplus, cezplus, dezplus;
			double aexbeyplus, bexaeyplus, bexceyplus, cexbeyplus;
			double cexdeyplus, dexceyplus, dexaeyplus, aexdeyplus;
			double aexceyplus, cexaeyplus, bexdeyplus, dexbeyplus;
			double det;
			double permanent, errbound;

			aex = pa[0] - pe[0];
			bex = pb[0] - pe[0];
			cex = pc[0] - pe[0];
			dex = pd[0] - pe[0];
			aey = pa[1] - pe[1];
			bey = pb[1] - pe[1];
			cey = pc[1] - pe[1];
			dey = pd[1] - pe[1];
			aez = pa[2] - pe[2];
			bez = pb[2] - pe[2];
			cez = pc[2] - pe[2];
			dez = pd[2] - pe[2];

			aexbey = aex * bey;
			bexaey = bex * aey;
			ab = aexbey - bexaey;
			bexcey = bex * cey;
			cexbey = cex * bey;
			bc = bexcey - cexbey;
			cexdey = cex * dey;
			dexcey = dex * cey;
			cd = cexdey - dexcey;
			dexaey = dex * aey;
			aexdey = aex * dey;
			da = dexaey - aexdey;

			aexcey = aex * cey;
			cexaey = cex * aey;
			ac = aexcey - cexaey;
			bexdey = bex * dey;
			dexbey = dex * bey;
			bd = bexdey - dexbey;

			abc = aez * bc - bez * ac + cez * ab;
			bcd = bez * cd - cez * bd + dez * bc;
			cda = cez * da + dez * ac + aez * cd;
			dab = dez * ab + aez * bd + bez * da;

			alift = aex * aex + aey * aey + aez * aez;
			blift = bex * bex + bey * bey + bez * bez;
			clift = cex * cex + cey * cey + cez * cez;
			dlift = dex * dex + dey * dey + dez * dez;

			det = (dlift * abc - clift * dab) + (blift * cda - alift * bcd);

			aezplus = Absolute(aez);
			bezplus = Absolute(bez);
			cezplus = Absolute(cez);
			dezplus = Absolute(dez);
			aexbeyplus = Absolute(aexbey);
			bexaeyplus = Absolute(bexaey);
			bexceyplus = Absolute(bexcey);
			cexbeyplus = Absolute(cexbey);
			cexdeyplus = Absolute(cexdey);
			dexceyplus = Absolute(dexcey);
			dexaeyplus = Absolute(dexaey);
			aexdeyplus = Absolute(aexdey);
			aexceyplus = Absolute(aexcey);
			cexaeyplus = Absolute(cexaey);
			bexdeyplus = Absolute(bexdey);
			dexbeyplus = Absolute(dexbey);
			permanent = ((cexdeyplus + dexceyplus) * bezplus
				+ (dexbeyplus + bexdeyplus) * cezplus
				+ (bexceyplus + cexbeyplus) * dezplus)
				* alift
				+ ((dexaeyplus + aexdeyplus) * cezplus
					+ (aexceyplus + cexaeyplus) * dezplus
					+ (cexdeyplus + dexceyplus) * aezplus)
				* blift
				+ ((aexbeyplus + bexaeyplus) * dezplus
					+ (bexdeyplus + dexbeyplus) * aezplus
					+ (dexaeyplus + aexdeyplus) * bezplus)
				* clift
				+ ((bexceyplus + cexbeyplus) * aezplus
					+ (cexaeyplus + aexceyplus) * bezplus
					+ (aexbeyplus + bexaeyplus) * cezplus)
				* dlift;
			errbound = isperrboundA * permanent;
			if ((det > errbound) || (-det > errbound)) {
				return det;
			}

			return insphereadapt(pa, pb, pc, pd, pe, permanent);
		}
	}

	double orient3d(std::array<Point3D, 4> const& points)
	{
		double adx, bdx, cdx, ady, bdy, cdy, adz, bdz, cdz;
		double bdxcdy, cdxbdy, cdxady, adxcdy, adxbdy, bdxady;
		double det;
		double permanent, errbound;

		adx = points[0].x - points[3].x;
		bdx = points[1].x - points[3].x;
		cdx = points[2].x - points[3].x;
		ady = points[0].y - points[3].y;
		bdy = points[1].y - points[3].y;
		cdy = points[2].y - points[3].y;
		adz = points[0].z - points[3].z;
		bdz = points[1].z - points[3].z;
		cdz = points[2].z - points[3].z;

		bdxcdy = bdx * cdy;
		cdxbdy = cdx * bdy;

		cdxady = cdx * ady;
		adxcdy = adx * cdy;

		adxbdy = adx * bdy;
		bdxady = bdx * ady;

		det = adz * (bdxcdy - cdxbdy)
			+ bdz * (cdxady - adxcdy)
			+ cdz * (adxbdy - bdxady);

		permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * Absolute(adz)
			+ (Absolute(cdxady) + Absolute(adxcdy)) * Absolute(bdz)
			+ (Absolute(adxbdy) + Absolute(bdxady)) * Absolute(cdz);
		errbound = o3derrboundA * permanent;
		if ((det > errbound) || (-det > errbound)) {
			return det;
		}

		return orient3dadapt(points,permanent);
	}

	/*double orient2d(boost::array<Point3D, 3> const& points)
	{
		double pa[2], pb[2], pc[2];

	}*/

	double insphere(std::array<Point3D, 5> const& points)
	{
		double pa[3], pb[3], pc[3], pd[3],pe[3];
		pa[0] = points[0].x;
		pa[1] = points[0].y;
		pa[2] = points[0].z;
		pb[0] = points[1].x;
		pb[1] = points[1].y;
		pb[2] = points[1].z;
		pc[0] = points[2].x;
		pc[1] = points[2].y;
		pc[2] = points[2].z;
		pd[0] = points[3].x;
		pd[1] = points[3].y;
		pd[2] = points[3].z;
		pe[0] = points[4].x;
		pe[1] = points[4].y;
		pe[2] = points[4].z;
		return insphere(pa, pb, pc, pd, pe);
	}
}
