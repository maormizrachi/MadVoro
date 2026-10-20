#ifndef MADVORO_PORTABILITY_MATH_HPP
#define MADVORO_PORTABILITY_MATH_HPP

#include <cmath>
#include <limits>

#include "KokkosTypes.hpp"

namespace MadVoro
{
namespace Portable
{

//! \brief std::numeric_limits<double>::min(), usable in device code.
constexpr double DOUBLE_MIN = std::numeric_limits<double>::min();

MADVORO_INLINE_FUNCTION double Sqrt(double x)
{
#ifdef MADVORO_WITH_GPU
    return Kokkos::sqrt(x);
#else
    return std::sqrt(x);
#endif
}

MADVORO_INLINE_FUNCTION double Abs(double x)
{
#ifdef MADVORO_WITH_GPU
    return Kokkos::abs(x);
#else
    return std::abs(x);
#endif
}

MADVORO_INLINE_FUNCTION double Pow(double x, double y)
{
#ifdef MADVORO_WITH_GPU
    return Kokkos::pow(x, y);
#else
    return std::pow(x, y);
#endif
}

//! \brief Same result as std::max(a, b), including for NaN arguments.
MADVORO_INLINE_FUNCTION double Max(double a, double b)
{
    return (a < b) ? b : a;
}

} // namespace Portable
} // namespace MadVoro

#endif // MADVORO_PORTABILITY_MATH_HPP
