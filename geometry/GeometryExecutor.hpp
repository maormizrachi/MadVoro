#ifndef MADVORO_GEOMETRY_EXECUTOR_HPP
#define MADVORO_GEOMETRY_EXECUTOR_HPP

// The executor of the geometry kernels: Kokkos (on a device or host threads) when
// MadVoro is built with Kokkos, the host one after the other otherwise. Both take the
// same calls and give the same results.
#ifdef MADVORO_WITH_KOKKOS
#include "KokkosGeometryExecutor.hpp"
#else
#include "HostGeometryExecutor.hpp"
#endif

namespace MadVoro
{
namespace Geometry
{
#ifdef MADVORO_WITH_KOKKOS
using GeometryExecutor = KokkosGeometryExecutor;
#else
using GeometryExecutor = HostGeometryExecutor;
#endif
} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_EXECUTOR_HPP
