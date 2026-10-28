#ifndef MADVORO_PORTABILITY_MATH_HPP
#define MADVORO_PORTABILITY_MATH_HPP

#include <cmath>
#include <limits>

#include "KokkosTypes.hpp"

namespace MadVoro
{
namespace Portable
{

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

// a * b + c with one rounding
MADVORO_INLINE_FUNCTION double Fma(double a, double b, double c)
{
#ifdef MADVORO_WITH_GPU
    return Kokkos::fma(a, b, c);
#else
    return std::fma(a, b, c);
#endif
}

// Same result as std::max(a, b), including for NaN arguments.
MADVORO_INLINE_FUNCTION double Max(double a, double b)
{
    return (a < b) ? b : a;
}

} // namespace Portable
} // namespace MadVoro

#endif // MADVORO_PORTABILITY_MATH_HPP
