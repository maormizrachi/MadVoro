#include <cstdint>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <mpi.h>

#include "regression_tests/lib/voronoi_test_common.hpp"

namespace
{

void Check(bool condition, int rank, int &failures, const std::string &message)
{
    if(condition)
    {
        return;
    }
    ++failures;
    std::cerr << "rank " << rank << ": " << message << std::endl;
}

int RunOwnershipChecks(MadVoro::regression_tests::VoronoiGrid &grid,
                       bool periodic,
                       std::uint64_t seed,
                       std::size_t numQueries,
                       int rank,
                       int nprocs)
{
    const Vector3D ll(0.0, 0.0, 0.0);
    const Vector3D ur(1.0, 1.0, 1.0);
    if(periodic)
    {
        grid.SetPeriodicBoundaries(true, true, true);
    }

    std::vector<Vector3D> meshPoints;
    if(rank == 0)
    {
        meshPoints = MadVoro::regression_tests::RandRectangular(
            nprocs > 1 ? static_cast<std::size_t>(20000) : static_cast<std::size_t>(5000),
            ll,
            ur,
            seed);
    }
    meshPoints = MadVoro::regression_tests::SpreadPointsFromRoot(meshPoints, rank);
    MadVoro::regression_tests::BuildVoronoiMesh(grid, meshPoints);

    int failures = 0;
    unsigned long long skippedQueries = 0;
    unsigned long long verifiedQueries = 0;
    const std::size_t localCells = grid.GetPointNo();
    for(std::size_t cell = 0; cell < localCells; ++cell)
    {
        const Vector3D center = grid.GetMeshPoint(cell);
        const std::size_t containing = grid.GetContainingCell(center);
        {
            std::ostringstream msg;
            msg << "owned center returned containing cell " << containing << " instead of " << cell;
            Check(containing == cell, rank, failures, msg.str());
        }
        Check(grid.IsPointInCell(center, cell), rank, failures, "owned center is not inside its local cell");
    }

    boost::mt19937_64 queryGen(seed + (periodic ? 1000003ULL : 17ULL));
    boost::random::uniform_real_distribution<> dist(-0.05, 1.05);
    for(std::size_t queryIdx = 0; queryIdx < numQueries; ++queryIdx)
    {
        Vector3D query(dist(queryGen), dist(queryGen), dist(queryGen));
        if(periodic)
        {
            grid.WrapPeriodicPoint(query);
        }
        else if(grid.IsPointOutsideBox(query))
        {
            ++skippedQueries;
            continue;
        }

        if(grid.GetOwner(query) != rank)
        {
            continue;
        }
        if(!grid.PointInMyDomain(query))
        {
            ++skippedQueries;
            continue;
        }

        const std::size_t containing = grid.GetContainingCell(query);
        const std::size_t bruteForce = MadVoro::regression_tests::FindOwnedCellBruteForce(grid, query);
        {
            std::ostringstream msg;
            msg << "query " << queryIdx << " GetContainingCell=" << containing
                << " brute_force=" << bruteForce;
            Check(bruteForce != static_cast<std::size_t>(-2), rank, failures, msg.str() + " (multiple owned matches)");
            Check(bruteForce == containing, rank, failures, msg.str());
        }
        Check(grid.IsPointInCell(query, containing), rank, failures,
              "GetContainingCell result failed IsPointInCell");
        ++verifiedQueries;
    }

    int globalFailures = 0;
    unsigned long long globalSkipped = 0;
    unsigned long long globalVerified = 0;
    MPI_Allreduce(&failures, &globalFailures, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&skippedQueries, &globalSkipped, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    MPI_Allreduce(&verifiedQueries, &globalVerified, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    if(rank == 0)
    {
        std::cout << (periodic ? "periodic" : "aperiodic")
                  << " verified_queries=" << globalVerified
                  << " skipped_queries=" << globalSkipped << std::endl;
    }

    return globalFailures;
}

} // namespace

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int rank = 0;
    int nprocs = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    const std::uint64_t seed = 987654321ULL;
    const std::size_t numQueries = 500;

    MadVoro::regression_tests::VoronoiGrid aperiodicGrid(Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 1.0, 1.0));
    const int aperiodicFailures = RunOwnershipChecks(aperiodicGrid, false, seed, numQueries, rank, nprocs);

    MadVoro::regression_tests::VoronoiGrid periodicGrid(Vector3D(0.0, 0.0, 0.0), Vector3D(1.0, 1.0, 1.0));
    const int periodicFailures = RunOwnershipChecks(periodicGrid, true, seed, numQueries, rank, nprocs);

    const int totalFailures = aperiodicFailures + periodicFailures;
    if(rank == 0)
    {
        if(totalFailures == 0)
        {
            std::cout << "voronoi_parallel_check PASS"
                      << " seed=" << seed
                      << " queries=" << numQueries
                      << " aperiodic_failures=" << aperiodicFailures
                      << " periodic_failures=" << periodicFailures
                      << std::endl;
        }
        else
        {
            std::cerr << "voronoi_parallel_check FAIL"
                      << " seed=" << seed
                      << " queries=" << numQueries
                      << " aperiodic_failures=" << aperiodicFailures
                      << " periodic_failures=" << periodicFailures
                      << std::endl;
        }
    }

    MPI_Finalize();
    return totalFailures == 0 ? 0 : 1;
}
