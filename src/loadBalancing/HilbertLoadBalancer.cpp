#ifdef MADVORO_WITH_MPI

#include "loadBalancing/HilbertLoadBalancer.hpp"
#include <iostream>

MadVoro::HilbertLoadBalancer::HilbertLoadBalancer(const std::shared_ptr<HilbertConvertor3D> convertor, const std::vector<hilbert_index_t> &boundaries)
    : CurveLoadBalancer(boundaries), convertor(convertor)
{}

void MadVoro::HilbertLoadBalancer::rebalance(const std::vector<Point3D> &points, const std::vector<double> &weights)
{
    if(this->convertor == nullptr)
    {
        throw MadVoro::Exception::MadVoroException("HilbertLoadBalancer::rebalance: convertor was not initialized yet");
    }

    std::vector<hilbert_index_t> indices;
    indices.reserve(points.size());
    for(const Point3D &point : points)
    {
        indices.push_back(this->convertor->xyz2d(point));
    }

    int dont_do_weights = (weights.empty() || std::all_of(weights.cbegin(), weights.cend(), [&weights](const double &x){return x == weights[0];}))? 1 : 0;
    MPI_Allreduce(MPI_IN_PLACE, &dont_do_weights, 1, MPI_INT, MPI_MAX, this->comm);
    if(this->rank == 0)
    {
        std::cout << "Running rebalancing" << std::endl;
    }
    if(dont_do_weights)
    {
        this->boundaries = MadVoro::MPI::Balance::getWeightedBorders2(indices, std::vector<double>(points.size(), 1.0));
    }
    else
    {
        this->boundaries = MadVoro::MPI::Balance::getWeightedBorders2(indices, weights);
    }
}

std::shared_ptr<MadVoro::LoadBalancer> MadVoro::HilbertLoadBalancer::clone(const std::shared_ptr<HilbertConvertor3D> newConvertor) const
{
    return std::make_shared<HilbertLoadBalancer>(newConvertor, this->boundaries);
}

#endif // MADVORO_WITH_MPI
