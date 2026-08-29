#ifndef MADVORO_VORONOI_TEST_COMMON_HPP
#define MADVORO_VORONOI_TEST_COMMON_HPP

#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_real_distribution.hpp>

#include "Voronoi3D.hpp"
#include "examples/Vector3D.hpp"

#ifdef MADVORO_WITH_MPI
#include <mpi.h>
#include <mpi_utils/mpi_collectives.hpp>
#endif

namespace MadVoro::regression_tests
{

using VoronoiGrid = MadVoro::Voronoi3D<Vector3D>;

inline std::vector<Vector3D> RandRectangular(std::size_t pointNum, const Vector3D &ll, const Vector3D &ur,
                                             boost::mt19937_64 &gen)
{
    boost::random::uniform_real_distribution<> dist;
    std::vector<Vector3D> points;
    points.reserve(pointNum);
    for(std::size_t i = 0; i < pointNum; ++i)
    {
        double x = ll.x + dist(gen) * (ur.x - ll.x);
        double y = ll.y + dist(gen) * (ur.y - ll.y);
        double z = ll.z + dist(gen) * (ur.z - ll.z);
        points.push_back(Vector3D(x, y, z));
    }
    return points;
}

inline std::vector<Vector3D> RandRectangular(std::size_t pointNum, const Vector3D &ll, const Vector3D &ur,
                                             std::uint64_t seed)
{
    boost::mt19937_64 gen(seed);
    return RandRectangular(pointNum, ll, ur, gen);
}

inline double SumOwnedCellVolumes(const VoronoiGrid &grid)
{
    double localVolume = 0.0;
    const std::size_t localCells = grid.GetPointNo();
    for(std::size_t cell = 0; cell < localCells; ++cell)
    {
        localVolume += grid.GetVolume(cell);
    }
    return localVolume;
}

inline double ReduceGlobalVolume(double localVolume)
{
#ifdef MADVORO_WITH_MPI
    double globalVolume = localVolume;
    MPI_Allreduce(&localVolume, &globalVolume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return globalVolume;
#else
    return localVolume;
#endif
}

inline double BoxVolume(const Vector3D &ll, const Vector3D &ur)
{
    return (ur.x - ll.x) * (ur.y - ll.y) * (ur.z - ll.z);
}

inline bool WriteVolumeMetrics(const std::string &filename, const std::string &modeLabel, bool periodic,
                               std::size_t numPoints, double totalVolume, double boxVolume, double relError,
                               int passFlag)
{
    std::ofstream out(filename);
    if(!out)
    {
        return false;
    }
    out.setf(std::ios::scientific);
    out.precision(16);
    out << "mode " << modeLabel << "\n";
    out << "periodic " << (periodic ? 1 : 0) << "\n";
    out << "num_points " << numPoints << "\n";
    out << "total_volume " << totalVolume << "\n";
    out << "box_volume " << boxVolume << "\n";
    out << "rel_error " << relError << "\n";
    out << "pass " << passFlag << "\n";
    return true;
}

#ifdef MADVORO_WITH_MPI
inline std::vector<Vector3D> SpreadPointsFromRoot(const std::vector<Vector3D> &rootPoints, int rank)
{
    std::vector<Vector3D> points;
    if(rank == 0)
    {
        points = rootPoints;
    }
    return MPI_Spread(points, 0, MPI_COMM_WORLD);
}
#endif

inline void BuildVoronoiMesh(VoronoiGrid &grid, const std::vector<Vector3D> &points)
{
#ifdef MADVORO_WITH_MPI
    grid.BuildParallel(points);
#else
    grid.Build(points);
#endif
}

inline std::size_t FindOwnedCellBruteForce(const VoronoiGrid &grid, const Vector3D &point)
{
    std::size_t foundCell = static_cast<std::size_t>(-1);
    const std::size_t localCells = grid.GetPointNo();
    for(std::size_t cell = 0; cell < localCells; ++cell)
    {
        if(!grid.IsPointInCell(point, cell))
        {
            continue;
        }
        if(foundCell != static_cast<std::size_t>(-1))
        {
            return static_cast<std::size_t>(-2);
        }
        foundCell = cell;
    }
    return foundCell;
}

} // namespace MadVoro::regression_tests

#endif // MADVORO_VORONOI_TEST_COMMON_HPP
