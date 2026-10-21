#ifndef MADVORO_GEOMETRY_EXECUTOR_HPP
#define MADVORO_GEOMETRY_EXECUTOR_HPP

// The executor of the geometry kernels: Kokkos (on a device or host threads) when
// MadVoro is built with Kokkos, the host one after the other otherwise. Both take the
// same calls and give the same results. MADVORO_GEOMETRY_ON_HOST keeps a Kokkos build
// on the host one, which then runs on Kokkos' host threads (the default execution
// space, e.g. a GPU, is left to other Kokkos users).
#if defined(MADVORO_WITH_GPU) && !defined(MADVORO_GEOMETRY_ON_HOST)
#include "KokkosGeometryExecutor.hpp"
#else
#include "HostGeometryExecutor.hpp"
#endif

namespace MadVoro
{
namespace Geometry
{
#if defined(MADVORO_WITH_GPU) && !defined(MADVORO_GEOMETRY_ON_HOST)
using GeometryExecutor = KokkosGeometryExecutor;
#else
using GeometryExecutor = HostGeometryExecutor;
#endif
} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_EXECUTOR_HPP
