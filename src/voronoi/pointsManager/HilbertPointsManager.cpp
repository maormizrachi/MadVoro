#ifdef MADVORO_WITH_MPI

#include "voronoi/pointsManager/HilbertPointsManager.hpp"

MadVoro::PointsExchangeResult MadVoro::HilbertPointsManager::exchange(const std::vector<Point3D> &allPoints, const std::vector<double> &allWeights, const std::vector<size_t> &indicesToWorkWith, const std::vector<double> &radiuses, const std::vector<Point3D> &previous_CM)
{
    PointsExchangeResult exchangeResult;
    if(this->envAgent != nullptr)
    {
        exchangeResult = this->pointsExchange([this](const _3DPointData &_point)
        {
            hilbert_index_t d = this->convertor->xyz2d(_point.point.x, _point.point.y, _point.point.z);
            size_t index = std::distance(this->responsibilityRange.cbegin(), std::upper_bound(this->responsibilityRange.cbegin(), this->responsibilityRange.cend(), d));
            return std::min<hilbert_index_t>(index, (this->size - 1));
        },
        allPoints, allWeights, indicesToWorkWith, radiuses, previous_CM); // exchange
        this->envAgent->updatePoints(exchangeResult.newPoints);
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

void MadVoro::HilbertPointsManager::rebalance(const std::vector<Point3D> &points, const std::vector<double> &weights)
{
    if(this->convertor == nullptr)
    {
        throw MadVoro::Exception::MadVoroException("HilbertPointsManager::rebalance: convertor was not initialized yet");
    }

    std::vector<hilbert_index_t> indices;
    for(const Point3D &point : points)
    {
        indices.push_back(this->convertor->xyz2d(point));
    }

    if(weights.empty())
    {
        this->responsibilityRange = MadVoro::MPI::Balance::getBorders(indices);
    }
    else
    {
        this->responsibilityRange = MadVoro::MPI::Balance::getWeightedBorders(indices, weights);
    }
    
    if(this->envAgent != nullptr)
    {
        this->envAgent->updateBorders(this->responsibilityRange, this->convertor->getOrder());
    }
}

/*
heuristic to determine the hilbert order.
Also initializes the convertor.
*/
void MadVoro::HilbertPointsManager::initializeHilbertParameters(const std::vector<Point3D> &points)
{
    DataStructure::OctTree<Point3D> tree(this->ll, this->ur, points);

    int depth = tree.getDepth(); // my own depth
    int hilbertOrder;
    MPI_Allreduce(&depth, &hilbertOrder, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD); // calculates maximal depth

    MPI_Allreduce(MPI_IN_PLACE, &this->ll.x, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &this->ll.y, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &this->ll.z, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &this->ur.x, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &this->ur.y, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE, &this->ur.z, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    // make a little bit space
    double x_length = this->ur.x - this->ll.x, y_length = this->ur.y - this->ll.y, z_length = this->ur.z - this->ll.z;
    this->ll.x -= std::abs(SPACE_FACTOR * x_length);
    this->ll.y -= std::abs(SPACE_FACTOR * y_length);
    this->ll.z -= std::abs(SPACE_FACTOR * z_length);
    this->ur.x += std::abs(SPACE_FACTOR * x_length);
    this->ur.y += std::abs(SPACE_FACTOR * y_length);
    this->ur.z += std::abs(SPACE_FACTOR * z_length);
    
    hilbertOrder = std::min<size_t>(MAX_HILBERT_ORDER, hilbertOrder);
    this->convertor = new HilbertRectangularConvertor3D(this->ll, this->ur, hilbertOrder);
}

MadVoro::PointsExchangeResult MadVoro::HilbertPointsManager::initialize(const std::vector<Point3D> &points, const std::vector<double> &weights, const std::vector<double> &radiuses, const std::vector<Point3D> &previous_CM)
{
    // if(this->rank == 0)
    // {
    //     std::cout << "initializes the points manager, and the environment agent" << std::endl;
    // }

    // calculate the first and initial order, and set it to the deepest hilbert order we have
    std::vector<size_t> allIndices(points.size());
    std::iota(allIndices.begin(), allIndices.end(), 0);

    this->initializeHilbertParameters(points); // also initializes the convertor

    this->rebalance(points, weights); // determines initial borders

    // making exchange according to these borders
    PointsExchangeResult exchangeResult = this->pointsExchange([this](const _3DPointData &_point)
    {
        hilbert_index_t d = this->convertor->xyz2d(_point.point.x, _point.point.y, _point.point.z);
        size_t index = std::distance(this->responsibilityRange.cbegin(), std::upper_bound(this->responsibilityRange.cbegin(), this->responsibilityRange.cend(), d));
        return std::min<hilbert_index_t>(index, (this->size - 1));
    },
    points, weights, allIndices, radiuses, previous_CM); // exchange
        
    this->envAgent = new HilbertTreeEnvironmentAgent(this->ll, this->ur, exchangeResult.newPoints, this->responsibilityRange, this->convertor);

    return exchangeResult;
}

#endif // MADVORO_WITH_MPI