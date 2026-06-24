#ifndef PREDICATES3D_HPP
#define PREDICATES3D_HPP 1

#include <vector>
#include "3D/elementary/Vector3D.hpp"
#include <array>

double orient3d(std::array<Vector3D,4> const& points);

double insphere(std::array<Vector3D, 5> const& points);

#endif //PREDICATES3D_HPP
