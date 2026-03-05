#ifndef CURVE_ENVIRONMENT_AGENT
#define CURVE_ENVIRONMENT_AGENT

#ifdef MADVORO_WITH_MPI

#include <vector>
#include <memory>
#include "loadBalancing/CurveLoadBalancer.hpp"
#include "EnvironmentAgent.h"

namespace MadVoro
{
    template<typename curve_index_t = hilbert_index_t, typename LoadBalancerType = CurveLoadBalancer>
    class CurveEnvironmentAgent : public EnvironmentAgent
    {
        static_assert(std::is_base_of<CurveLoadBalancer, LoadBalancerType>::value, "LoadBalancerType must inherit from CurveLoadBalancer");

    public:
        inline CurveEnvironmentAgent(const Point3D &ll, const Point3D &ur, const std::shared_ptr<LoadBalancerType> curveLoadBalancer, const MPI_Comm &comm = MPI_COMM_WORLD):
            EnvironmentAgent(ll, ur, comm), loadBalancer(curveLoadBalancer)
        {}

        virtual ~CurveEnvironmentAgent() = default;

        virtual void setLoadBalancer(std::shared_ptr<LoadBalancerType> newLoadBalancer)
        {
            this->loadBalancer = newLoadBalancer;
            this->onRebalance();
        }

        virtual inline int getCellOwner(curve_index_t d) const
        {
            size_t index = static_cast<size_t>(std::distance(this->loadBalancer->boundaries.cbegin(), std::upper_bound(this->loadBalancer->boundaries.cbegin(), this->loadBalancer->boundaries.cend(), d)));
            return std::min<int>(static_cast<int>(index), this->size - 1);
        };

        virtual void onExchange(const std::vector<Point3D> &newPoints) override
        {}

        virtual void onRebalance(void) override
        {}

    protected:
        std::shared_ptr<LoadBalancerType> loadBalancer;
    };
}

#endif // MADVORO_WITH_MPI

#endif // CURVE_ENVIRONMENT_AGENT
