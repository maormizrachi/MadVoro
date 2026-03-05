#ifndef POINTS_MANAGER_HPP
#define POINTS_MANAGER_HPP

#ifdef MADVORO_WITH_MPI

#include <algorithm>
#include <numeric>
#include <vector>
#include <cmath>
#include <memory>
#include <mpi.h>
#include <assert.h>

#include "mpi/balance/balance.hpp"
#include "mpi/exchange/exchange.hpp"
#include "elementary/Point3D.hpp"
#include "environment/EnvironmentAgent.h"
#include "loadBalancing/LoadBalancer.hpp"

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
        std::vector<bool> participatingIndices;
    };

    struct PointData : public MadVoro::MPI::Serializable
    {
        size_t indexInAllPoints;
        Point3D point;
        double radius;
        Point3D CM;
        double weight;
        bool participating;

        PointData() = default;

        size_t dump(MadVoro::MPI::Serializer *serializer) const override
        {
            size_t bytes = 0;
            bytes += serializer->insert(this->indexInAllPoints);
            bytes += serializer->insert(this->point);
            bytes += serializer->insert(this->radius);
            bytes += serializer->insert(this->CM);
            bytes += serializer->insert(this->weight);
            bytes += serializer->insert(this->participating);
            return bytes;
        }

        size_t load(const MadVoro::MPI::Serializer *serializer, size_t byteOffset) override
        {
            size_t bytes = 0;
            bytes += serializer->extract(this->indexInAllPoints, byteOffset);
            bytes += serializer->extract(this->point, byteOffset + bytes);
            bytes += serializer->extract(this->radius, byteOffset + bytes);
            bytes += serializer->extract(this->CM, byteOffset + bytes);
            bytes += serializer->extract(this->weight, byteOffset + bytes);
            bytes += serializer->extract(this->participating, byteOffset + bytes);
            return bytes;
        }
    };

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

        virtual std::shared_ptr<PointsManager> clone(void) const = 0;

        PointsManager &operator=(const PointsManager &other) = delete;

        virtual PointsExchangeResult exchange(const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToWorkWith, const std::vector<double> &radiuses, const std::vector<Point3D> &previous_CM, bool noExchange = false) = 0;

        virtual void rebalance(const std::vector<Point3D> &points, const std::vector<double> &weights = std::vector<double>()) = 0;

        virtual const std::shared_ptr<EnvironmentAgent> getEnvironmentAgent() const = 0;

        virtual void setLoadBalancer(std::shared_ptr<LoadBalancer> loadBalancer) = 0;

        virtual std::shared_ptr<LoadBalancer> getLoadBalancer(void) = 0;

        void setImbalanceTolerance(double tolerance);

        void reportImbalance(void) const;

        bool checkForRebalance(double myWeight) const;

        bool shouldRebalance(const std::vector<double> &weights) const;

        inline bool shouldRebalance(void) const{return this->checkForRebalance(this->totalWeight);};

        PointsExchangeResult update(const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToWorkWith, const std::vector<double> &radiuses, const std::vector<Point3D> &previous_CM, bool doRebalance = true, bool doExchange = true);

    protected:
        Point3D ll, ur;
        MPI_Comm comm;
        int rank, size;
        double totalWeight;
        double imbalanceTolerance = BALANCE_FACTOR;

        /**
         * performs a point exchange, according to a given determination function (point -> rank)
        */
        template<typename DetermineFunc>
        PointsExchangeResult pointsExchange(const DetermineFunc &func, const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToWorkWith, const std::vector<double> &radiuses, const std::vector<Point3D> &previous_CM) const;
    };

    template<typename DetermineFunc>
    PointsExchangeResult PointsManager::pointsExchange(const DetermineFunc &func, const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToWorkWith, const std::vector<double> &radiuses, const std::vector<Point3D> &previous_CM) const
    {
        std::vector<PointData> data;
        data.reserve(allPoints.size());
        for(size_t pointIdx = 0; pointIdx < allPoints.size(); pointIdx++)
        {
            const Point3D &point = allPoints[pointIdx];
            data.emplace_back();
            PointData &pointRadius = data.back();
            pointRadius.indexInAllPoints = pointIdx;
            pointRadius.point = point;
            pointRadius.radius = radiuses[pointIdx];
            pointRadius.weight = allWeights[pointIdx];
            pointRadius.CM = previous_CM[pointIdx];
            pointRadius.participating = false;
        }
        
        for(const size_t &pointIdx : indicesToWorkWith)
        {
            data[pointIdx].participating = true;
        }

        ExchangeAnswer<PointData> answer = MadVoro::MPI::Exchange::dataExchange(data, func, this->comm);

        PointsExchangeResult toReturn;
        
        toReturn.indicesToSelf = std::move(answer.indicesToMe);
        toReturn.sentProcessors = std::move(answer.processesSend);
        toReturn.sentIndicesToProcessors = std::move(answer.indicesToProcesses);

        std::vector<PointData> &ans = answer.output;
        toReturn.newPoints.reserve(ans.size());
        toReturn.newRadiuses.reserve(ans.size());
        toReturn.newCMs.reserve(ans.size());
        toReturn.newWeights.reserve(ans.size());
        toReturn.participatingIndices.resize(ans.size(), false);

        for(const PointData &_point : ans)
        {
            size_t pointIdx = toReturn.newPoints.size();
            toReturn.newPoints.push_back(_point.point);
            if(_point.participating)
            {
                toReturn.participatingIndices[pointIdx] = true;
            }
            toReturn.newRadiuses.push_back(_point.radius);
            toReturn.newCMs.push_back(_point.CM);
            toReturn.newWeights.push_back(_point.weight);
        }

        assert(toReturn.newPoints.size() == toReturn.newWeights.size());
        return toReturn;
    };
}

#endif // MADVORO_WITH_MPI

#endif // POINTS_MANAGER_HPP
