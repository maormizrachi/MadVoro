#ifndef _HILBERT_TREE_ENVIRONMENT_AGENT_HPP
#define _HILBERT_TREE_ENVIRONMENT_AGENT_HPP

#ifdef MADVORO_WITH_MPI

#include "HilbertCurveEnvAgent.hpp"
#include "hilbert/rectangular/HilbertRectangularTree3D.hpp"

namespace MadVoro
{
    class HilbertTreeEnvironmentAgent : public HilbertCurveEnvironmentAgent
    {
    public:
        using HilbertTree_Type = MadVoro::DataStructure::HilbertTree3D<DEFAULT_RANKS_IN_LEAVES>;

        inline HilbertTreeEnvironmentAgent(const Point3D &ll, const Point3D &ur, const std::shared_ptr<HilbertLoadBalancer> loadBalancer, const MPI_Comm &comm = MPI_COMM_WORLD): 
                HilbertCurveEnvironmentAgent(ll, ur, loadBalancer, comm)
        {
            this->rectangularConvertor = dynamic_cast<const HilbertRectangularConvertor3D*>(this->loadBalancer->convertor.get());
            if(this->rectangularConvertor == nullptr)
            {
                throw MadVoro::Exception::MadVoroException("'HilbertTreeEnvironmentAgent' should be initialized with a rectangular hilbert convertor");
            }
            this->hilbertTree = std::make_shared<HilbertTree_Type>(this->rectangularConvertor, this->loadBalancer->boundaries, this->comm);
        };

        inline ~HilbertTreeEnvironmentAgent() override = default;

        inline std::shared_ptr<HilbertCurveEnvironmentAgent> clone(const std::shared_ptr<HilbertLoadBalancer> newLoadBalancer) const override
        {
            return std::make_shared<HilbertTreeEnvironmentAgent>(this->ll, this->ur, newLoadBalancer, this->comm);
        }

        inline EnvironmentAgent::RanksSet getIntersectingRanks(const Point3D &center, double radius) const override
        {
            return this->hilbertTree->getIntersectingRanks(center, radius);
        };

        inline void onExchange(const std::vector<Point3D> &newPoints) override
        {
            this->HilbertCurveEnvironmentAgent::onExchange(newPoints);
        }
        
        inline void onRebalance(void) override
        {
            this->HilbertCurveEnvironmentAgent::onRebalance();
            this->hilbertTree = std::make_shared<HilbertTree_Type>(
                dynamic_cast<const HilbertRectangularConvertor3D*>(this->loadBalancer->convertor.get()),
                this->loadBalancer->boundaries, this->comm);
        }
        
        template<typename U>
        inline HilbertCurveEnvironmentAgent::DistancesVector getClosestFurthestPointsByRanks(const U &point) const
        {
            return this->hilbertTree->getClosestFurthestPointsByRanks(point);
        }

        inline const std::shared_ptr<HilbertTree_Type> &getHilbertTree() const{return this->hilbertTree;};

    private:
        std::shared_ptr<HilbertTree_Type> hilbertTree;
        const HilbertRectangularConvertor3D *rectangularConvertor;
    };
}

#endif // MADVORO_WITH_MPI

#endif // _HILBERT_TREE_ENVIRONMENT_AGENT_HPP
