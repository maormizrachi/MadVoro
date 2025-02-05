#ifndef HILBERT_CURVE_ENVAGENT_HPP
#define HILBERT_CURVE_ENVAGENT_HPP

#ifdef MADVORO_WITH_MPI

#include "../CurveEnvAgent.hpp"
#include "hilbert/HilbertConvertor3D.hpp"

namespace MadVoro
{
    class HilbertCurveEnvironmentAgent : public CurveEnvironmentAgent<hilbert_index_t>
    {
    public:
        using DistancesVector = std::vector<std::pair<typename Point3D::coord_type, typename Point3D::coord_type>>;

        inline HilbertCurveEnvironmentAgent(const Point3D &ll, const Point3D &ur, const std::vector<hilbert_index_t> &ranges, HilbertConvertor3D *convertor, const MPI_Comm &comm = MPI_COMM_WORLD):
            CurveEnvironmentAgent(ll, ur, ranges, comm), convertor(convertor)
        {
            this->order = this->convertor->getOrder();
        };

        virtual ~HilbertCurveEnvironmentAgent() = default;

        virtual inline int getOwner(const Point3D &point) const override
        {
            return this->getCellOwner(this->convertor->xyz2d(point));
        };

        virtual void updatePoints(const std::vector<Point3D> &newPoints) override
        {}

        virtual inline void updateBorders(const std::vector<hilbert_index_t> &newRange, int newOrder)
        {
            this->CurveEnvironmentAgent::updateBorders(newRange);
            if(this->convertor != nullptr)
            {
                this->convertor->changeOrder(newOrder);
            }
        }

    protected:
        HilbertConvertor3D *convertor = nullptr;
        int order;
    };
}

#endif // MADVORO_WITH_MPI

#endif // HILBERT_CURVE_ENVAGENT_HPP
