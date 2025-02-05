#ifndef PREDICATES3D_HPP
#define PREDICATES3D_HPP 1

#include <vector>
#include "elementary/Point3D.hpp"
#include <array>

namespace MadVoro
{
    #if defined(__clang__) || defined (__GNUC__)
    # define ATTRIBUTE_NO_SANITIZE_ADDRESS __attribute__((no_sanitize_address))
    #else
    # define ATTRIBUTE_NO_SANITIZE_ADDRESS
    #endif
    ATTRIBUTE_NO_SANITIZE_ADDRESS
    double orient3d(std::array<Point3D,4> const& points);

    double insphere(std::array<Point3D, 5> const& points);
}

#endif //PREDICATES3D_HPP
