#ifndef PREDICATES3D_HPP
#define PREDICATES3D_HPP 1

#include <array>
#include <cmath>
#include "Predicates3D_internal.hpp"

namespace MadVoro {

template <typename PointT>
double orient3d(std::array<PointT, 4> const& points)
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

	permanent = (std::abs(bdxcdy) + std::abs(cdxbdy)) * std::abs(adz)
		+ (std::abs(cdxady) + std::abs(adxcdy)) * std::abs(bdz)
		+ (std::abs(adxbdy) + std::abs(bdxady)) * std::abs(cdz);
	errbound = o3derrboundA * permanent;
	if ((det > errbound) || (-det > errbound)) {
		return det;
	}

	return orient3dadapt(points, permanent);
}

template <typename PointT>
double insphere(std::array<PointT, 5> const& points)
{
	double pa[3], pb[3], pc[3], pd[3], pe[3];
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

} // namespace MadVoro

#endif // PREDICATES3D_HPP
