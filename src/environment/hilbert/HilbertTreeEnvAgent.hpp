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

        inline HilbertTreeEnvironmentAgent(const Point3D &ll, const Point3D &ur, const std::vector<Point3D> &points, const std::vector<hilbert_index_t> &ranges, HilbertConvertor3D *convertor, const MPI_Comm &comm = MPI_COMM_WORLD): 
                HilbertCurveEnvironmentAgent(ll, ur, ranges, convertor, comm)
        {
            this->rectangularConvertor = dynamic_cast<HilbertRectangularConvertor3D*>(this->convertor);
            if(this->rectangularConvertor == nullptr)
            {
                throw MadVoro::Exception::MadVoroException("'HilbertTreeEnvironmentAgent' should be initialized with a rectangular hilbert convertor");
            }

            this->hilbertTree = new HilbertTree_Type(this->rectangularConvertor, this->range, this->comm);
        };

        inline ~HilbertTreeEnvironmentAgent() override
        {
            delete this->hilbertTree;
        };

        inline EnvironmentAgent::RanksSet getIntersectingRanks(const Point3D &center, double radius) const override
        {
            return this->hilbertTree->getIntersectingRanks(center, radius);
        };

        inline void updatePoints(const std::vector<Point3D> &newPoints) override
        {
            this->HilbertCurveEnvironmentAgent::updatePoints(newPoints);
        }
        
        inline void updateBorders(const std::vector<hilbert_index_t> &newRange, int newOrder) override
        {
            this->HilbertCurveEnvironmentAgent::updateBorders(newRange, newOrder);
            delete this->hilbertTree;
            this->hilbertTree = new HilbertTree_Type(this->rectangularConvertor, this->range, this->comm);
        }

        inline int getOrder() const{return this->order;};
        
        template<typename U>
        inline HilbertCurveEnvironmentAgent::DistancesVector getClosestFurthestPointsByRanks(const U &point) const
        {
            return this->hilbertTree->getClosestFurthestPointsByRanks(point);
        }

        inline const HilbertTree_Type *getHilbertTree() const{return this->hilbertTree;};

    private:
        const HilbertTree_Type *hilbertTree;
        HilbertRectangularConvertor3D *rectangularConvertor;
    };
}

#endif // MADVORO_WITH_MPI

#endif // _HILBERT_TREE_ENVIRONMENT_AGENT_HPP