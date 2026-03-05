#ifdef MADVORO_WITH_MPI

#include "loadBalancing/CurveLoadBalancer.hpp"

MadVoro::CurveLoadBalancer::CurveLoadBalancer(const std::vector<hilbert_index_t> &boundaries)
    : LoadBalancer(), boundaries(boundaries)
{}

#endif // MADVORO_WITH_MPI
