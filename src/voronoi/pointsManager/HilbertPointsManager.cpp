#ifdef MADVORO_WITH_MPI

#include "voronoi/pointsManager/HilbertPointsManager.hpp"
#include <memory>

MadVoro::HilbertPointsManager::HilbertPointsManager(const Point3D &ll, const Point3D &ur, const MPI_Comm &comm)
    : PointsManager(ll, ur, comm)
{}

std::shared_ptr<MadVoro::PointsManager> MadVoro::HilbertPointsManager::clone(void) const
{
    auto cloned = std::make_shared<HilbertPointsManager>(this->ll, this->ur, this->comm);

    if(this->convertor != nullptr)
    {
        cloned->convertor = this->convertor->clone();
    }

    if(this->loadBalancer != nullptr)
    {
        cloned->loadBalancer = std::dynamic_pointer_cast<HilbertLoadBalancer>(
            this->loadBalancer->clone(cloned->convertor));
    }

    if(this->envAgent != nullptr)
    {
        cloned->envAgent = this->envAgent->clone(cloned->loadBalancer);
    }
    return cloned;
}

MadVoro::PointsExchangeResult MadVoro::HilbertPointsManager::exchange(const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToWorkWith, const std::vector<double> &radiuses, const std::vector<Point3D> &previous_CM, bool noExchange)
{
    PointsExchangeResult exchangeResult;

    if(this->envAgent != nullptr)
    {
        const std::vector<hilbert_index_t> &responsibilityRange = this->loadBalancer->boundaries;

        if(noExchange)
        {
            exchangeResult = this->pointsExchange([this](const PointData &_point)
            {
                return this->rank;
            },
            allPoints, allWeights, indicesToWorkWith, radiuses, previous_CM);
        }
        else
        {
            exchangeResult = this->pointsExchange([this, &responsibilityRange](const PointData &_point)
            {
                hilbert_index_t d = this->convertor->xyz2d(_point.point.x, _point.point.y, _point.point.z);
                size_t index = std::distance(responsibilityRange.cbegin(), std::upper_bound(responsibilityRange.cbegin(), responsibilityRange.cend(), d));
                return std::min<hilbert_index_t>(index, (this->size - 1));
            },
            allPoints, allWeights, indicesToWorkWith, radiuses, previous_CM);
        }
        this->envAgent->onExchange(exchangeResult.newPoints);
    }
    else
    {
        if(allPoints.size() != indicesToWorkWith.size())
        {
            MadVoro::Exception::MadVoroException eo("Error in HilbertPointsManager::exchange: in the first build, a mesh with all the points should be built. Currently, points and indicesToWorkWith have different sizes");
            eo.addEntry("allPoints.size()", allPoints.size());
            eo.addEntry("indicesToWorkWith.size()", indicesToWorkWith.size());
            throw eo;
        }
        exchangeResult = this->initialize(allPoints, allWeights, radiuses, previous_CM);
    }

    return exchangeResult;
}

void MadVoro::HilbertPointsManager::setLoadBalancer(std::shared_ptr<LoadBalancer> newLoadBalancer)
{
    HilbertLoadBalancer *hilbertLoadBalancer = dynamic_cast<HilbertLoadBalancer*>(newLoadBalancer.get());
    if(hilbertLoadBalancer == nullptr)
    {
        throw MadVoro::Exception::MadVoroException("HilbertPointsManager::setLoadBalancer: given load balancer is not a HilbertLoadBalancer");
    }
    if(this->rank == 0)
    {
        std::cout << "Restoring Load Balancer" << std::endl;
    }

    this->loadBalancer = std::dynamic_pointer_cast<HilbertLoadBalancer>(newLoadBalancer);
    this->loadBalancer->convertor = this->convertor;
    this->envAgent->setLoadBalancer(this->loadBalancer);
}

std::shared_ptr<MadVoro::LoadBalancer> MadVoro::HilbertPointsManager::getLoadBalancer(void)
{
    return this->loadBalancer;
}

void MadVoro::HilbertPointsManager::rebalance(const std::vector<Point3D> &points, const std::vector<double> &weights)
{
    this->loadBalancer = std::dynamic_pointer_cast<HilbertLoadBalancer>(
        this->loadBalancer->clone(this->convertor));
    this->loadBalancer->rebalance(points, weights);
    if(this->envAgent != nullptr)
    {
        this->envAgent->setLoadBalancer(this->loadBalancer);
    }
}

void MadVoro::HilbertPointsManager::initializeHilbertParameters(const std::vector<Point3D> &points)
{
    DataStructure::OctTree<Point3D> tree(this->ll, this->ur, points);

    int depth = tree.getDepth();
    int hilbertOrder;
    MPI_Allreduce(&depth, &hilbertOrder, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

    MPI_Allreduce(MPI_IN_PLACE, &this->ll.x, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &this->ll.y, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &this->ll.z, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &this->ur.x, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &this->ur.y, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &this->ur.z, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    double x_length = this->ur.x - this->ll.x, y_length = this->ur.y - this->ll.y, z_length = this->ur.z - this->ll.z;
    this->ll.x -= std::abs(SPACE_FACTOR * x_length);
    this->ll.y -= std::abs(SPACE_FACTOR * y_length);
    this->ll.z -= std::abs(SPACE_FACTOR * z_length);
    this->ur.x += std::abs(SPACE_FACTOR * x_length);
    this->ur.y += std::abs(SPACE_FACTOR * y_length);
    this->ur.z += std::abs(SPACE_FACTOR * z_length);
    
    hilbertOrder = std::min<size_t>(MAX_HILBERT_ORDER, hilbertOrder);
    this->convertor = std::make_shared<HilbertRectangularConvertor3D>(this->ll, this->ur, hilbertOrder);
}

MadVoro::PointsExchangeResult MadVoro::HilbertPointsManager::initialize(const std::vector<Point3D> &points, const std::vector<double> &weights, const std::vector<double> &radiuses, const std::vector<Point3D> &previous_CM)
{
    std::vector<size_t> allIndices(points.size());
    std::iota(allIndices.begin(), allIndices.end(), 0);

    this->initializeHilbertParameters(points);

    this->loadBalancer = std::make_shared<HilbertLoadBalancer>(this->convertor);

    this->rebalance(points, weights);

    const std::vector<hilbert_index_t> &responsibilityRange = this->loadBalancer->boundaries;

    PointsExchangeResult exchangeResult = this->pointsExchange([this, &responsibilityRange](const PointData &_point)
    {
        hilbert_index_t d = this->convertor->xyz2d(_point.point.x, _point.point.y, _point.point.z);
        size_t index = std::distance(responsibilityRange.cbegin(), std::upper_bound(responsibilityRange.cbegin(), responsibilityRange.cend(), d));
        return std::min<size_t>(index, (this->size - 1));
    },
    points, weights, allIndices, radiuses, previous_CM);
        
    this->envAgent = std::make_shared<HilbertTreeEnvironmentAgent>(this->ll, this->ur, this->loadBalancer, this->comm);

    return exchangeResult;
}

#endif // MADVORO_WITH_MPI
