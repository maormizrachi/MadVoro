#ifdef MADVORO_WITH_MPI

#include "PointsManager.hpp"

using namespace MadVoro;

bool MadVoro::PointsManager::checkForRebalance(double myWeight) const
{
    // checks if I have too many weight, and notify other ranks
    double totalWeight;
    MPI_Allreduce(&myWeight, &totalWeight, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    double idealWeight = totalWeight / this->size;
    int I_say = (myWeight >= (BALANCE_FACTOR * idealWeight))? 1 : 0; // if I say 'rebalance' or not
    int rebalance = 0; // if someone says 'rebalance' or not
    MPI_Allreduce(&I_say, &rebalance, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
    if((rebalance > 0) and (this->rank == 0))
    {
        std::cout << "doing rebalance" << std::endl;
    }
    return (rebalance > 0);
}

MadVoro::PointsExchangeResult MadVoro::PointsManager::update(const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToWorkWith, const std::vector<double> &radiuses, const std::vector<Point3D> &previous_CM, bool allowRebalance)
{
    // if envAgent is null, the `exchange` will perform an initialization as well.
    // `rebalance` is used only when the environment agent is initialized.
    PointsExchangeResult result = this->exchange(allPoints, allWeights, indicesToWorkWith, radiuses, previous_CM);
    this->totalWeight = std::accumulate(result.newWeights.cbegin(), result.newWeights.cend(), 0.0);
    if(allowRebalance and this->checkForRebalance(this->totalWeight))
    {
        assert(this->getEnvironmentAgent() != nullptr);
        this->rebalance(allPoints, allWeights);
        result = this->exchange(allPoints, allWeights, indicesToWorkWith, radiuses, previous_CM);
        this->totalWeight = std::accumulate(result.newWeights.cbegin(), result.newWeights.cend(), 0.0);
    }
    // std::cout << "total weight of rank " << this->rank << " is " << this->totalWeight << " with " << result.newPoints.size() << " points" << std::endl;
    return result;
}

#endif // MADVORO_WITH_MPI