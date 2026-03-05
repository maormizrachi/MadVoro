#ifndef CURVE_LOAD_BALANCER_HPP
#define CURVE_LOAD_BALANCER_HPP

#include "hilbert/hilbertTypes.h"

#ifdef MADVORO_WITH_MPI

#include <vector>
#include "LoadBalancer.hpp"

namespace MadVoro
{
    class CurveLoadBalancer : public LoadBalancer
    {
    public:
        CurveLoadBalancer(const std::vector<hilbert_index_t> &boundaries = std::vector<hilbert_index_t>());

        virtual ~CurveLoadBalancer() override = default;

        std::vector<hilbert_index_t> boundaries;
    };
}

#endif // MADVORO_WITH_MPI

#endif // CURVE_LOAD_BALANCER_HPP
