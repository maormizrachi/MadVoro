#ifndef _HILBERT_ENVIRONMENT_AGENT_HPP
#define _HILBERT_ENVIRONMENT_AGENT_HPP

#ifdef MADVORO_WITH_MPI

#include "HilbertCurveEnvAgent.hpp"

#define AVERAGE_INTERSECT 128
#define NULL_ORDER -1

namespace MadVoro
{
    class HilbertEnvironmentAgent : public HilbertCurveEnvironmentAgent
    {
    public:
        using CellsSet = boost::container::flat_set<hilbert_index_t>;

        inline HilbertEnvironmentAgent(const Point3D &ll, const Point3D &ur, const std::shared_ptr<HilbertLoadBalancer> loadBalancer, const MPI_Comm &comm = MPI_COMM_WORLD):
                HilbertCurveEnvironmentAgent(ll, ur, loadBalancer, comm)
        {
            this->setOrder(this->getOrder());
        };
            
        RanksSet getIntersectingRanks(const Point3D &center, double radius) const override;
        
        CellsSet getIntersectingCells(const Point3D &center, double radius) const;

        inline std::shared_ptr<HilbertCurveEnvironmentAgent> clone(const std::shared_ptr<HilbertLoadBalancer> newLoadBalancer) const override
        {
            return std::make_shared<HilbertEnvironmentAgent>(this->ll, this->ur, newLoadBalancer, this->comm);
        }

        inline void onExchange(const std::vector<Point3D> &newPoints) override
        {
            this->HilbertCurveEnvironmentAgent::onExchange(newPoints);
        }

        inline void onRebalance(void) override
        {
            this->HilbertCurveEnvironmentAgent::onRebalance();
            this->setOrder(this->getOrder());
        }

    private:
        Point3D sidesLengths;

        inline void setOrder(int order)
        {
            if(order == NULL_ORDER)
            {
                return;
            }
            int clampedOrder = std::min<int>(order, MAX_HILBERT_ORDER);
            this->sidesLengths = (this->ur - this->ll) / pow(2, clampedOrder);
            this->loadBalancer->convertor->changeOrder(clampedOrder);
        }
    };
}

#endif // MADVORO_WITH_MPI

#endif // _HILBERT_ENVIRONMENT_AGENT_HPP
