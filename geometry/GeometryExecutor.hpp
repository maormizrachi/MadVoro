#ifndef MADVORO_GEOMETRY_EXECUTOR_HPP
#define MADVORO_GEOMETRY_EXECUTOR_HPP

// MADVORO_GEOMETRY_ON_HOST keeps a Kokkos build on the host executor (on Kokkos' host threads).
#if defined(MADVORO_WITH_GPU) and not defined(MADVORO_GEOMETRY_ON_HOST)
#include "KokkosGeometryExecutor.hpp"
#else
#include "HostGeometryExecutor.hpp"
#endif

namespace MadVoro
{
namespace Geometry
{
#if defined(MADVORO_WITH_GPU) and not defined(MADVORO_GEOMETRY_ON_HOST)
using GeometryExecutor = KokkosGeometryExecutor;
#else
using GeometryExecutor = HostGeometryExecutor;
#endif
} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_EXECUTOR_HPP
