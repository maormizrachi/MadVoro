#ifndef INTERSECTIONS_HPP
#define INTERSECTIONS_HPP 1

#include "ds/utils/Sphere.hpp"
#include "elementary/Face3D.hpp"

namespace MadVoro
{
    bool Face3DSphereIntersections(Face3D const& face, MadVoro::Geometry::Sphere<Point3D> const& sphere, Point3D const& normal);
}

#endif //INTERSECTIONS_HPP
