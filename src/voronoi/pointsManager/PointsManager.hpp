#ifndef POINTS_MANAGER_HPP
#define POINTS_MANAGER_HPP

#ifdef MADVORO_WITH_MPI

#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>
#include <mpi.h>
#include <assert.h>

#include "mpi/balance/balance.hpp"
#include "mpi/exchange/exchange.hpp"
#include "utils/point/3DPoint.hpp"
#include "elementary/Point3D.hpp"
#include "environment/EnvironmentAgent.h"

#define BALANCE_FACTOR 1.15

namespace MadVoro
{
    /**
     * \author Maor Mizrachi
     * \brief A result for points exchange running.
    */
    struct PointsExchangeResult
    {
        std::vector<Point3D> newPoints;
        std::vector<double> newRadiuses;
        std::vector<int> sentProcessors;
        std::vector<std::vector<size_t>> sentIndicesToProcessors;
        std::vector<size_t> indicesToSelf;
        std::vector<Point3D> newCMs;
        std::vector<double> newWeights;
        std::vector<size_t> participatingIndices;
    };

    typedef struct _3DPointData
    {
        size_t indexInAllPoints;
        _3DPoint point;
        double radius;
        _3DPoint CM;
        double weight;
        bool participating;
    } _3DPointData;

    /**
     * \author Maor Mizrachi
     * \brief A point manager performs data movement between ranks (borders determination and points exchange according to borders).
    */
    class PointsManager
    {
    public:
        inline PointsManager(const Point3D &ll, const Point3D &ur, const MPI_Comm &comm = MPI_COMM_WORLD): ll(ll), ur(ur), comm(comm), totalWeight(0)
        {
            MPI_Comm_size(this->comm, &this->size);
            MPI_Comm_rank(this->comm, &this->rank);
        };

        virtual ~PointsManager() = default;

        PointsManager &operator=(const PointsManager &other) = delete;

        virtual PointsExchangeResult exchange(const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToWorkWith, const std::vector<double> &radiuses, const std::vector<Point3D> &previous_CM) = 0;

        virtual void rebalance(const std::vector<Point3D> &points, const std::vector<double> &weights = std::vector<double>()) = 0;

        virtual const EnvironmentAgent *getEnvironmentAgent() const = 0;

        bool checkForRebalance(double myWeight) const;

        PointsExchangeResult update(const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToWorkWith, const std::vector<double> &radiuses, const std::vector<Point3D> &previous_CM, bool allowRebalance = true);

    protected:
        Point3D ll, ur;
        MPI_Comm comm;
        int rank, size;
        double totalWeight;

        /**
         * performs a point exchange, according to a given determination function (point -> rank)
        */
        template<typename DetermineFunc>
        PointsExchangeResult pointsExchange(const DetermineFunc &func, const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToWorkWith, const std::vector<double> &radiuses, const std::vector<Point3D> &previous_CM) const;
    };

    template<typename DetermineFunc>
    PointsExchangeResult PointsManager::pointsExchange(const DetermineFunc &func, const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToWorkWith, const std::vector<double> &radiuses, const std::vector<Point3D> &previous_CM) const
    {
        std::vector<_3DPointData> data;
        data.reserve(allPoints.size());
        for(size_t pointIdx = 0; pointIdx < allPoints.size(); pointIdx++)
        {
            const Point3D &point = allPoints[pointIdx];
            data.emplace_back();
            _3DPointData &pointRadius = data.back();
            pointRadius.indexInAllPoints = pointIdx;
            pointRadius.point = _3DPoint(point.x, point.y, point.z);
            pointRadius.radius = radiuses[pointIdx];
            pointRadius.weight = allWeights[pointIdx];
            pointRadius.CM = _3DPoint(previous_CM[pointIdx].x, previous_CM[pointIdx].y, previous_CM[pointIdx].z);
            pointRadius.participating = false;
        }
        
        for(const size_t &pointIdx : indicesToWorkWith)
        {
            data[pointIdx].participating = true;
        }

        // // re-build the function so that it maintains the points that are not participating
        // auto new_func = [&func, this, &participating](const _3DPointData &point){return ((not participating[point.indexInAllPoints])? this->rank : func(point));};
        ExchangeAnswer<_3DPointData> answer = MadVoro::MPI::Exchange::dataExchange(data, func, this->comm);

        // arrange the return value data structure
        PointsExchangeResult toReturn;
        
        toReturn.indicesToSelf = std::move(answer.indicesToMe);
        toReturn.sentProcessors = std::move(answer.processesSend);
        toReturn.sentIndicesToProcessors = std::move(answer.indicesToProcesses);

        std::vector<_3DPointData> &ans = answer.output;
        toReturn.newPoints.reserve(ans.size());
        toReturn.newRadiuses.reserve(ans.size());
        toReturn.newCMs.reserve(ans.size());
        toReturn.newWeights.reserve(ans.size());
        toReturn.participatingIndices.resize(ans.size(), false);

        for(const _3DPointData &_point : ans)
        {
            size_t pointIdx = toReturn.newPoints.size();
            toReturn.newPoints.emplace_back(_point.point.x, _point.point.y, _point.point.z);
            if(_point.participating)
            {
                toReturn.participatingIndices[pointIdx] = true;
            }
            toReturn.newRadiuses.push_back(_point.radius);
            toReturn.newCMs.emplace_back(_point.CM.x, _point.CM.y, _point.CM.z);
            toReturn.newWeights.push_back(_point.weight);
        }

        assert(toReturn.newPoints.size() == toReturn.newWeights.size());
        return toReturn;
    };
}

#endif // MADVORO_WITH_MPI

#endif // POINTS_MANAGER_HPP