#ifndef MADVORO_KOKKOS_TYPES_HPP
#define MADVORO_KOKKOS_TYPES_HPP

#ifdef MADVORO_WITH_KOKKOS
#include <Kokkos_Core.hpp>

#define MADVORO_FUNCTION KOKKOS_FUNCTION
#define MADVORO_INLINE_FUNCTION KOKKOS_INLINE_FUNCTION
#else
#define MADVORO_FUNCTION
#define MADVORO_INLINE_FUNCTION inline
#endif

// Before a host-device function template that some callers instantiate with host-only
// types (e.g. a std::vector sink): tells nvcc not to warn about those host instantiations.
#if defined(__CUDACC__)
#define MADVORO_HOST_DEVICE_TEMPLATE _Pragma("nv_exec_check_disable")
#else
#define MADVORO_HOST_DEVICE_TEMPLATE
#endif

#endif // MADVORO_KOKKOS_TYPES_HPP
