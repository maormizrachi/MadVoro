#include <cmath>
#include <cstddef>
#include <iostream>
#include <string>

#include "regression_tests/lib/voronoi_test_common.hpp"

#ifdef MADVORO_WITH_MPI
#include <mpi.h>
#endif

int main()
{
    int rank = 0;
    int worldSize = 1;
#ifdef MADVORO_WITH_MPI
    MPI_Init(nullptr, nullptr);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
#endif

    const Vector3D ll(0.0, 0.0, 0.0);
    const Vector3D ur(1.0, 1.0, 1.0);
    const std::uint64_t seed = 424242ULL;
    const std::size_t numPoints = (worldSize > 1)
        ? static_cast<std::size_t>(100000)
        : static_cast<std::size_t>(10000);

    std::vector<Vector3D> points;
#ifdef MADVORO_WITH_MPI
    if(rank == 0)
    {
        points = MadVoro::regression_tests::RandRectangular(numPoints, ll, ur, seed);
    }
    points = MadVoro::regression_tests::SpreadPointsFromRoot(points, rank);
#else
    points = MadVoro::regression_tests::RandRectangular(numPoints, ll, ur, seed);
#endif

    MadVoro::regression_tests::VoronoiGrid grid(ll, ur);
    MadVoro::regression_tests::BuildVoronoiMesh(grid, points);

    const double totalVolume = MadVoro::regression_tests::ReduceGlobalVolume(
        MadVoro::regression_tests::SumOwnedCellVolumes(grid));
    const double boxVolume = MadVoro::regression_tests::BoxVolume(ll, ur);
    const double relError = std::abs(totalVolume - boxVolume) / boxVolume;
    const int passed = (relError < 1e-10) ? 1 : 0;

    if(rank == 0)
    {
        std::cout << "voronoi_volume seed=" << seed << "\n"
                  << "total_volume = " << totalVolume << "\n"
                  << "box_volume   = " << boxVolume << "\n"
                  << "rel_error    = " << relError << "\n"
                  << "pass         = " << passed << std::endl;

        if(!MadVoro::regression_tests::WriteVolumeMetrics("voronoi_volume_metrics.txt",
                                                          worldSize > 1 ? "mpi" : "serial",
                                                          false,
                                                          numPoints,
                                                          totalVolume,
                                                          boxVolume,
                                                          relError,
                                                          passed))
        {
            std::cerr << "voronoi_volume FAIL: could not write metrics file" << std::endl;
#ifdef MADVORO_WITH_MPI
            MPI_Abort(MPI_COMM_WORLD, 2);
#endif
            return 2;
        }
    }

#ifdef MADVORO_WITH_MPI
    MPI_Finalize();
#endif
    return passed ? 0 : 1;
}
