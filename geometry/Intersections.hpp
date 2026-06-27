#ifndef INTERSECTIONS_HPP
#define INTERSECTIONS_HPP 1

#include <spatial_ds/utils/Sphere.hpp>
#include "../elementary/Face3D.hpp"
#include "../elementary/PointOps.hpp"

namespace MadVoro {

using namespace MadVoro::fallback;

namespace {

template <typename PointT>
bool PointInPolygon(Face3D<PointT> const& face, PointT const& point)
{
	PointT normal = CrossProduct(face.vertices[0] - point, face.vertices[1] - point);
	const size_t Nloop = face.vertices.size() - 1;
	for (size_t i = 0; i < Nloop; ++i)
		if (ScalarProd(CrossProduct(face.vertices[i + 1] - point, face.vertices[(i + 2) % (Nloop + 1)] - point),
			normal) < 0)
			return false;
	return true;
}

template <typename PointT>
bool CircleSegmentIntersect(PointT const& p0, PointT const& p1, PointT const& center, double R)
{
	PointT AC = center - p0;
	PointT AB = p1 - p0;
	double d = ScalarProd(AC, AB);
	if (d < 0)
	{
		if (fastabs(AC) > R)
			return false;
		else
			return true;
	}
	double LAB = fastabs(AB);
	if (d > LAB * LAB)
	{
		if (fastabs(center - p1) > R)
			return false;
		else
			return true;
	}
	PointT closest = p0 + AB * d / (LAB * LAB);
	if (fastabs(center - closest) > R)
		return false;
	else
		return true;
}

} // anonymous namespace

template <typename PointT>
bool FaceSphereIntersections(Face3D<PointT> const& face, Sphere<PointT> const& sphere, PointT const& normal)
{
	double D = ScalarProd(normal, sphere.center - face.vertices[0]);

	if (std::abs(D) > sphere.radius)
		return false;
	PointT circle_center;
	circle_center.x = sphere.center.x - D * normal.x;
	circle_center.y = sphere.center.y - D * normal.y;
	circle_center.z = sphere.center.z - D * normal.z;
	std::size_t Nloop = face.vertices.size();
	if (PointInPolygon(face, circle_center))
		return true;
	double R = std::sqrt(sphere.radius * sphere.radius - D * D);
	for (std::size_t i = 0; i < Nloop; ++i)
	{
		if (CircleSegmentIntersect(face.vertices[(i + 1) % Nloop], face.vertices[i], circle_center, R))
			return true;
	}
	return false;
}

} // namespace MadVoro

#endif // INTERSECTIONS_HPP
