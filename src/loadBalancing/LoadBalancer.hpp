#ifndef LOAD_BALANCER_HPP
#define LOAD_BALANCER_HPP

#ifdef MADVORO_WITH_MPI

#include <vector>
#include <mpi.h>
#include "elementary/Point3D.hpp"
#include "mpi/serialize/mpi_commands.hpp"

namespace MadVoro
{
    class LoadBalancer
    {
    public:
        LoadBalancer(const MPI_Comm &comm = MPI_COMM_WORLD);

        virtual void rebalance(const std::vector<Point3D> &points, const std::vector<double> &weights = std::vector<double>()) = 0;

        virtual ~LoadBalancer() = default;

    protected:
        MPI_Comm comm;
        rank_t rank, size;
    };
}

#endif // MADVORO_WITH_MPI

#endif // LOAD_BALANCER_HPP
