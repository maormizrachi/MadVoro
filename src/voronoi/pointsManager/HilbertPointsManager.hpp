#ifndef HILBERT_POINTS_MANAGER_HPP
#define HILBERT_POINTS_MANAGER_HPP

#ifdef MADVORO_WITH_MPI

#include <numeric>
#include <memory>
#include "PointsManager.hpp"
#include "environment/hilbert/DistributedOctEnvAgent.hpp"
#include "environment/hilbert/HilbertTreeEnvAgent.hpp"
#include "hilbert/rectangular/HilbertRectangularConvertor3D.hpp"
#include "hilbert/ordinary/HilbertOrdinaryConvertor3D.hpp"

#define SPACE_FACTOR 1e-5

namespace MadVoro
{
    class HilbertPointsManager : public PointsManager
    {
    public:
        HilbertPointsManager(const Point3D &ll, const Point3D &ur, const MPI_Comm &comm = MPI_COMM_WORLD)
            : PointsManager(ll, ur, comm), envAgent(nullptr), convertor(nullptr)
        {}

        inline ~HilbertPointsManager() override
        {
            delete this->envAgent;
            delete this->convertor;
        };

        inline const EnvironmentAgent *getEnvironmentAgent() const override
        {
            return this->envAgent;
        }

        HilbertPointsManager &operator=(const HilbertPointsManager &other) = delete;

        PointsExchangeResult exchange(const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToWorkWith, const std::vector<double> &radiuses, const std::vector<Point3D> &previous_CM) override;

        void rebalance(const std::vector<Point3D> &points, const std::vector<double> &weights = std::vector<double>()) override;
        
    private:
        void initializeHilbertParameters(const std::vector<Point3D> &points);

        PointsExchangeResult initialize(const std::vector<Point3D> &points, const std::vector<double> &weights, const std::vector<double> &radiuses, const std::vector<Point3D> &previous_CM);

        HilbertCurveEnvironmentAgent *envAgent;
        HilbertConvertor3D *convertor;
        std::vector<hilbert_index_t> responsibilityRange;
    };
}

#endif // MADVORO_WITH_MPI

#endif // HILBERT_POINTS_MANAGER_HPP