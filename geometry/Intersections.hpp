#ifndef INTERSECTIONS_HPP
#define INTERSECTIONS_HPP 1

#include <spatial_ds/utils/Sphere.hpp>
#include "3D/elementary/Face.hpp"

bool FaceSphereIntersections(Face const& face, Sphere<Vector3D> const& sphere, Vector3D const& normal);

#endif //INTERSECTIONS_HPP
