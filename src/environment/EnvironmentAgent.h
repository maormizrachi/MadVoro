#ifndef ENVIRONMENT_MADVORO_H
#define ENVIRONMENT_MADVORO_H

#ifdef MADVORO_WITH_MPI

#include <mpi.h>
#include <boost/container/flat_set.hpp>

namespace MadVoro
{
    /**
     * \author Maor Mizrachi
     * \brief The environment agent is responsible for "knowing" the environment. It can calculate the ranks that intersect a sphere, or calculate the owner of a certain point.
    */
    class EnvironmentAgent
    {
    public:
        using RanksSet = boost::container::flat_set<int>;
        
        inline EnvironmentAgent(const Point3D &ll, const Point3D &ur, const MPI_Comm &comm = MPI_COMM_WORLD): ll(ll), ur(ur), comm(comm)
        {
            MPI_Comm_rank(this->comm, &this->rank);
            MPI_Comm_size(this->comm, &this->size);
        };

        virtual ~EnvironmentAgent() = default;
        
        virtual RanksSet getIntersectingRanks(const Point3D &center, double radius) const = 0;

        virtual int getOwner(const Point3D &point) const = 0;

    protected:
        Point3D ll, ur;
        MPI_Comm comm;
        int rank, size;
    };
}

#endif // MADVORO_WITH_MPI

#endif // ENVIRONMENT_MADVORO_H