#include <cmath>
#include <cstddef>
#include <iostream>

#include "regression_tests/lib/voronoi_test_common.hpp"

#ifdef MADVORO_WITH_MPI
#include <mpi.h>
#endif

namespace
{

bool CellHasPhysicalBoundary(const MadVoro::regression_tests::VoronoiGrid &grid, std::size_t cellIndex)
{
    const MadVoro::face_vec &faces = grid.GetCellFaces(cellIndex);
    for(std::size_t faceIdx : faces)
    {
        if(grid.BoundaryFace(faceIdx))
        {
            return true;
        }
    }
    return false;
}

bool ApproxEqual(double a, double b, double relTol = 1e-6)
{
    double scale = std::max(std::abs(a), std::abs(b));
    if(scale < 1.0)
    {
        scale = 1.0;
    }
    return std::abs(a - b) <= relTol * scale;
}

bool TestOnePointPeriodicSerial(MadVoro::regression_tests::VoronoiGrid &grid)
{
    grid.SetPeriodicBoundaries(true, true, true);
    std::vector<Vector3D> points = {Vector3D(0.5, 0.5, 0.5)};
    grid.Build(points);

    const double boxVolume = 1.0;
    const double cellVolume = grid.GetVolume(0);
    if(!ApproxEqual(cellVolume, boxVolume))
    {
        std::cerr << "One-point periodic: expected cell volume " << boxVolume
                  << ", got " << cellVolume << std::endl;
        return false;
    }
    if(CellHasPhysicalBoundary(grid, 0))
    {
        std::cerr << "One-point periodic: cell has physical boundary faces" << std::endl;
        return false;
    }
    return true;
}

bool TestLatticePeriodicSerial(MadVoro::regression_tests::VoronoiGrid &grid)
{
    grid.SetPeriodicBoundaries(true, true, true);
    std::vector<Vector3D> points;
    for(int ix = 0; ix < 2; ++ix)
    {
        for(int iy = 0; iy < 2; ++iy)
        {
            for(int iz = 0; iz < 2; ++iz)
            {
                points.emplace_back(0.25 + 0.5 * ix, 0.25 + 0.5 * iy, 0.25 + 0.5 * iz);
            }
        }
    }
    grid.Build(points);

    const double expectedCellVolume = 1.0 / 8.0;
    for(std::size_t i = 0; i < points.size(); ++i)
    {
        if(!ApproxEqual(grid.GetVolume(i), expectedCellVolume))
        {
            std::cerr << "Lattice periodic: cell " << i << " volume mismatch" << std::endl;
            return false;
        }
        if(CellHasPhysicalBoundary(grid, i))
        {
            std::cerr << "Lattice periodic: cell " << i << " has physical boundary faces" << std::endl;
            return false;
        }
    }
    return true;
}

} // namespace

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
    const std::uint64_t seed = 515151ULL;
    const std::size_t numPoints = (worldSize > 1)
        ? static_cast<std::size_t>(100000)
        : static_cast<std::size_t>(10000);

    int passed = 1;
#ifndef MADVORO_WITH_MPI
    {
        MadVoro::regression_tests::VoronoiGrid serialGrid(ll, ur);
        passed = TestOnePointPeriodicSerial(serialGrid) ? passed : 0;
        MadVoro::regression_tests::VoronoiGrid latticeGrid(ll, ur);
        passed = TestLatticePeriodicSerial(latticeGrid) ? passed : 0;
    }
#endif

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
    grid.SetPeriodicBoundaries(true, true, true);
    MadVoro::regression_tests::BuildVoronoiMesh(grid, points);

    const double totalVolume = MadVoro::regression_tests::ReduceGlobalVolume(
        MadVoro::regression_tests::SumOwnedCellVolumes(grid));
    const double boxVolume = MadVoro::regression_tests::BoxVolume(ll, ur);
    const double relError = std::abs(totalVolume - boxVolume) / boxVolume;
    passed = (relError < 1e-10) ? passed : 0;

    if(rank == 0)
    {
        std::cout << "voronoi_volume_periodic seed=" << seed << "\n"
                  << "total_volume = " << totalVolume << "\n"
                  << "box_volume   = " << boxVolume << "\n"
                  << "rel_error    = " << relError << "\n"
                  << "pass         = " << passed << std::endl;

        if(!MadVoro::regression_tests::WriteVolumeMetrics("voronoi_volume_periodic_metrics.txt",
                                                          worldSize > 1 ? "mpi" : "serial",
                                                          true,
                                                          numPoints,
                                                          totalVolume,
                                                          boxVolume,
                                                          relError,
                                                          passed))
        {
            std::cerr << "voronoi_volume_periodic FAIL: could not write metrics file" << std::endl;
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
