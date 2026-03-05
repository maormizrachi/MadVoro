#ifdef MADVORO_WITH_MPI

#include "PointsManager.hpp"
#include <chrono>

using namespace MadVoro;

void MadVoro::PointsManager::setImbalanceTolerance(double tolerance)
{
    this->imbalanceTolerance = tolerance;
}

void MadVoro::PointsManager::reportImbalance(void) const
{
    struct
    {
        double weight;
        int rank;
    } myWeightedRank, maxWeighted, minWeighted;
    myWeightedRank.weight = this->totalWeight;
    myWeightedRank.rank = this->rank;
    MPI_Allreduce(&myWeightedRank, &maxWeighted, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    MPI_Allreduce(&myWeightedRank, &minWeighted, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);
    double avgWeight;
    MPI_Allreduce(&this->totalWeight, &avgWeight, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    avgWeight /= static_cast<double>(this->size);
    if(this->rank == 0)
    {
        std::cout << "Imbalance report: max weight is " << maxWeighted.weight << " in rank " << maxWeighted.rank << ", min weight is " << minWeighted.weight << " in rank " << minWeighted.rank << ", average weight is " << avgWeight << std::endl;
    }
}

bool MadVoro::PointsManager::checkForRebalance(double myWeight) const
{
    struct
    {
        double weight;
        int rank;
    } myWeightRanked, maxWeight;

    double totalWeight;
    MPI_Allreduce(&myWeight, &totalWeight, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double idealWeight = totalWeight / this->size;
    myWeightRanked.weight = myWeight;
    myWeightRanked.rank = this->rank;
    MPI_Allreduce(&myWeightRanked, &maxWeight, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
    if(this->rank == 0)
    {
        std::cout << "Max weight in rank " << maxWeight.rank << ", its weight is " << maxWeight.weight << ", ideal is " << idealWeight << std::endl;
    }
    if(maxWeight.weight >= (this->imbalanceTolerance * idealWeight))
    {
        if(this->rank == 0)
        {
            std::cout << "Doing rebalance!" << std::endl;
        }
        return true;
    }
    return false;
}

bool MadVoro::PointsManager::shouldRebalance(const std::vector<double> &weights) const
{
    double totalWeight = std::accumulate(weights.cbegin(), weights.cend(), 0.0);
    return this->checkForRebalance(totalWeight);
}

MadVoro::PointsExchangeResult MadVoro::PointsManager::update(const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToWorkWith, const std::vector<double> &radiuses, const std::vector<Point3D> &previous_CM, bool doRebalance, bool doExchange)
{
    PointsExchangeResult result = this->exchange(allPoints, allWeights, indicesToWorkWith, radiuses, previous_CM, not doExchange);
    this->totalWeight = std::accumulate(result.newWeights.cbegin(), result.newWeights.cend(), 0.0);
    if(doRebalance and this->checkForRebalance(this->totalWeight))
    {
        assert(this->getEnvironmentAgent() != nullptr);
        this->rebalance(allPoints, allWeights);
        result = this->exchange(allPoints, allWeights, indicesToWorkWith, radiuses, previous_CM, not doExchange);
        this->totalWeight = std::accumulate(result.newWeights.cbegin(), result.newWeights.cend(), 0.0);
        this->reportImbalance();
    }
    return result;
}

#endif // MADVORO_WITH_MPI
