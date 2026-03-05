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
#include "loadBalancing/HilbertLoadBalancer.hpp"

#define SPACE_FACTOR 1e-5

namespace MadVoro
{
    class HilbertPointsManager : public PointsManager
    {
    public:
        HilbertPointsManager(const Point3D &ll, const Point3D &ur, const MPI_Comm &comm = MPI_COMM_WORLD);

        inline ~HilbertPointsManager() override = default;

        std::shared_ptr<PointsManager> clone(void) const override;

        inline const std::shared_ptr<EnvironmentAgent> getEnvironmentAgent() const override{return this->envAgent;}

        HilbertPointsManager &operator=(const HilbertPointsManager &other) = delete;

        PointsExchangeResult exchange(const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToWorkWith, const std::vector<double> &radiuses, const std::vector<Point3D> &previous_CM, bool noExchange = false) override;

        void rebalance(const std::vector<Point3D> &points, const std::vector<double> &weights = std::vector<double>()) override;
        
        void setLoadBalancer(std::shared_ptr<LoadBalancer> loadBalancer) override;

        std::shared_ptr<LoadBalancer> getLoadBalancer(void) override;

    private:
        void initializeHilbertParameters(const std::vector<Point3D> &points);

        PointsExchangeResult initialize(const std::vector<Point3D> &points, const std::vector<double> &weights, const std::vector<double> &radiuses, const std::vector<Point3D> &previous_CM);

        std::shared_ptr<HilbertLoadBalancer> loadBalancer = nullptr;
        std::shared_ptr<HilbertCurveEnvironmentAgent> envAgent = nullptr;
        std::shared_ptr<HilbertConvertor3D> convertor = nullptr;
    };
}

#endif // MADVORO_WITH_MPI

#endif // HILBERT_POINTS_MANAGER_HPP
