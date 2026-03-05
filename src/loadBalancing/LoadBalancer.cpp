#ifdef MADVORO_WITH_MPI

#include "loadBalancing/LoadBalancer.hpp"

MadVoro::LoadBalancer::LoadBalancer(const MPI_Comm &comm)
    : comm(comm)
{
    MPI_Comm_rank(this->comm, &this->rank);
    MPI_Comm_size(this->comm, &this->size);
}

#endif // MADVORO_WITH_MPI
