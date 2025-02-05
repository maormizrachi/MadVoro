#ifndef CURVE_ENVIRONMENT_AGENT
#define CURVE_ENVIRONMENT_AGENT

#ifdef MADVORO_WITH_MPI

#include <vector>
#include "EnvironmentAgent.h"

namespace MadVoro
{
    template<typename curve_index_t = size_t>
    class CurveEnvironmentAgent : public EnvironmentAgent
    {
    public:
        inline CurveEnvironmentAgent(const Point3D &ll, const Point3D &ur, const std::vector<curve_index_t> &ranges, const MPI_Comm &comm = MPI_COMM_WORLD):
            EnvironmentAgent(ll, ur, comm), range(ranges)
        {}

        virtual ~CurveEnvironmentAgent() = default;

        virtual inline int getCellOwner(curve_index_t d) const
        {
            int index = static_cast<int>(std::distance(this->range.begin(), std::upper_bound(this->range.begin(), this->range.end(), d)));
            return std::min<int>(index, this->size - 1);
        };

        virtual void updatePoints(const std::vector<Point3D> &newPoints)
        {}

        virtual inline void updateBorders(const std::vector<curve_index_t> &newRange)
        {
            this->range = newRange;
        }

    protected:
        std::vector<curve_index_t> range;
    };
}

#endif // MADVORO_WITH_MPI

#endif // CURVE_ENVIRONMENT_AGENT