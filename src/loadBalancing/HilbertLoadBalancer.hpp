#ifndef HILBERT_LOAD_BALANCER_HPP
#define HILBERT_LOAD_BALANCER_HPP

#ifdef MADVORO_WITH_MPI

#include "CurveLoadBalancer.hpp"
#include "hilbert/HilbertConvertor3D.hpp"
#include "mpi/balance/balance.hpp"
#include "mpi/balance/weightedBalance2.hpp"
#include <memory>

namespace MadVoro
{
    class HilbertLoadBalancer : public CurveLoadBalancer
    {
    public:
        HilbertLoadBalancer(const std::shared_ptr<HilbertConvertor3D> convertor, const std::vector<hilbert_index_t> &boundaries = std::vector<hilbert_index_t>());

        ~HilbertLoadBalancer() override = default;

        void rebalance(const std::vector<Point3D> &points, const std::vector<double> &weights) override;

        std::shared_ptr<LoadBalancer> clone(const std::shared_ptr<HilbertConvertor3D> newConvertor) const;

        std::shared_ptr<HilbertConvertor3D> convertor;
    };
}

#endif // MADVORO_WITH_MPI

#endif // HILBERT_LOAD_BALANCER_HPP
