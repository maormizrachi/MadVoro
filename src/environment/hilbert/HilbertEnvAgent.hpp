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

        inline HilbertEnvironmentAgent(const Point3D &ll, const Point3D &ur, const std::vector<hilbert_index_t> &ranges, HilbertConvertor3D *convertor, const MPI_Comm &comm = MPI_COMM_WORLD):
                HilbertCurveEnvironmentAgent(ll, ur, ranges, convertor, comm)
        {
            this->setOrder(this->order);
        };
            
        RanksSet getIntersectingRanks(const Point3D &center, double radius) const override;
        
        CellsSet getIntersectingCells(const Point3D &center, double radius) const;
        
        inline int getOrder() const{return this->order;};

        inline void updatePoints(const std::vector<Point3D> &newPoints) override
        {
            this->HilbertCurveEnvironmentAgent::updatePoints(newPoints);
        }

        inline void updateBorders(const std::vector<hilbert_index_t> &newRange, int newOrder) override
        {
            this->range = newRange;
            this->setOrder(newOrder);
        }

    private:
        Point3D sidesLengths;
        HilbertConvertor3D *convertor = nullptr;

        inline void setOrder(int order)
        {
            if(order == NULL_ORDER)
            {
                return;
            }
            this->order = std::min<int>(order, MAX_HILBERT_ORDER);
            this->sidesLengths = (this->ur - this->ll) / pow(2, this->order);
            this->convertor->changeOrder(this->order);
        }
    };
}

#endif // MADVORO_WITH_MPI

#endif // _HILBERT_ENVIRONMENT_AGENT_HPP