#ifndef _DIST_OCT_ENVIRONMENT_AGENT_HPP
#define _DIST_OCT_ENVIRONMENT_AGENT_HPP

#ifdef MADVORO_WITH_MPI

#include "HilbertCurveEnvAgent.hpp"
#include "ds/DistributedOctTree/DistributedOctTree.hpp"

#define RANKS_IN_LEAF 4

namespace MadVoro
{
    class DistributedOctEnvironmentAgent : public HilbertCurveEnvironmentAgent
    {
    public:
        using DistributedOctTree_Type = MadVoro::DataStructure::DistributedOctTree<Point3D, RANKS_IN_LEAF>;

        inline DistributedOctEnvironmentAgent(const Point3D &ll, const Point3D &ur, const std::vector<Point3D> &points, const std::vector<hilbert_index_t> &ranges, const MPI_Comm &comm = MPI_COMM_WORLD): 
                HilbertCurveEnvironmentAgent(ll, ur, ranges, convertor, comm)
        {
            DataStructure::OctTree<Point3D> myTree(this->ll, this->ur, points);
            this->distributedOctTree = new DistributedOctTree_Type(&myTree, false /* no detailed nodes info */, this->comm);
        };

        inline ~DistributedOctEnvironmentAgent(){delete this->distributedOctTree;};

        inline EnvironmentAgent::RanksSet getIntersectingRanks(const Point3D &center, double radius) const override
        {
            return this->distributedOctTree->getIntersectingRanks(center, radius);
        };

        inline void updatePoints(const std::vector<Point3D> &newPoints) override
        {
            this->HilbertCurveEnvironmentAgent::updatePoints(newPoints);
            delete this->distributedOctTree;
            DataStructure::OctTree<Point3D> myTree(this->ll, this->ur, newPoints);
            this->distributedOctTree = new DistributedOctTree_Type(&myTree, false, this->comm);
        }

        inline int getOwner(const Point3D &point) const override
        {
            return this->getCellOwner(this->convertor->xyz2d(point));
        };

        inline void updateBorders(const std::vector<hilbert_index_t> &newRange, int newOrder) override
        {
            this->HilbertCurveEnvironmentAgent::updateBorders(newRange, newOrder);
        }

        const DistributedOctTree_Type *getOctTree() const{return this->distributedOctTree;};

        inline int getOrder() const{return this->order;};
        
        template<typename U>
        inline HilbertCurveEnvironmentAgent::DistancesVector getClosestFurthestPointsByRanks(const U &point) const
        {
            return this->distributedOctTree->getClosestFurthestPointsByRanks(point);
        }

    private:
        DistributedOctTree_Type *distributedOctTree = nullptr;
        HilbertConvertor3D *convertor = nullptr;
        std::vector<hilbert_index_t> range;
        int order;
    };
}

#endif // MADVORO_WITH_MPI

#endif // _DIST_OCT_ENVIRONMENT_AGENT_HPP