#ifndef HILBERT_CURVE_ENVAGENT_HPP
#define HILBERT_CURVE_ENVAGENT_HPP

#ifdef MADVORO_WITH_MPI

#include "../CurveEnvAgent.hpp"
#include "hilbert/HilbertConvertor3D.hpp"
#include "loadBalancing/HilbertLoadBalancer.hpp"

namespace MadVoro
{
    class HilbertCurveEnvironmentAgent : public CurveEnvironmentAgent<hilbert_index_t, HilbertLoadBalancer>
    {
    public:
        using DistancesVector = std::vector<std::pair<typename Point3D::coord_type, typename Point3D::coord_type>>;

        inline HilbertCurveEnvironmentAgent(const Point3D &ll, const Point3D &ur, const std::shared_ptr<HilbertLoadBalancer> loadBalancer, const MPI_Comm &comm = MPI_COMM_WORLD):
            CurveEnvironmentAgent<hilbert_index_t, HilbertLoadBalancer>(ll, ur, loadBalancer, comm)
        {};

        virtual ~HilbertCurveEnvironmentAgent() = default;

        virtual std::shared_ptr<HilbertCurveEnvironmentAgent> clone(const std::shared_ptr<HilbertLoadBalancer> newLoadBalancer) const = 0;

        virtual inline int getOwner(const Point3D &point) const override
        {
            return this->getCellOwner(this->loadBalancer->convertor->xyz2d(point));
        };

        virtual void onExchange(const std::vector<Point3D> &newPoints) override
        {
            this->CurveEnvironmentAgent::onExchange(newPoints);
        }

        virtual void onRebalance(void) override
        {
            this->CurveEnvironmentAgent::onRebalance();
        }

        inline int getOrder() const{return static_cast<int>(this->loadBalancer->convertor->getOrder());};
    };
}

#endif // MADVORO_WITH_MPI

#endif // HILBERT_CURVE_ENVAGENT_HPP
