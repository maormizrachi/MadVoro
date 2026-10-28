// Building the cells of a Cartesian grid by clipping: every cell is the cube around its point (volume h^3, centroid at the point, faces with its 6 axis neighbors), as the Delaunay builder's, in a first build and in rebuilds of the same points, in a box and periodic
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "regression_tests/lib/voronoi_test_common.hpp"

#ifdef MADVORO_WITH_MPI
#include <mpi.h>
#endif

namespace
{

using MadVoro::regression_tests::VoronoiGrid;
using Key = std::array<double, 3>;

// relative to the cell's size: the cubes' moments to round-off, and faces this small do not count
constexpr double MOMENTS_TOLERANCE = 1e-10;
constexpr double TINY_FACE = 1e-8;
constexpr int BUILDS = 4;

Key KeyOf(const Vector3D &p)
{
    return {p.x, p.y, p.z};
}

std::size_t GlobalSum(std::size_t value)
{
#ifdef MADVORO_WITH_MPI
    unsigned long long sum = value;
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    return static_cast<std::size_t>(sum);
#else
    return value;
#endif
}

double GlobalMax(double value)
{
#ifdef MADVORO_WITH_MPI
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
#endif
    return value;
}

// the centers of side^3 cubes filling the unit box (rank 0's; spread over the ranks)
std::vector<Vector3D> MakeGrid(std::size_t side, int rank)
{
    std::vector<Vector3D> points;
    if(rank == 0)
    {
        const double h = 1.0 / static_cast<double>(side);
        for(std::size_t i = 0; i < side; ++i)
        {
            for(std::size_t j = 0; j < side; ++j)
            {
                for(std::size_t k = 0; k < side; ++k)
                {
                    points.emplace_back(h * (i + 0.5), h * (j + 0.5), h * (k + 0.5));
                }
            }
        }
    }
    return MadVoro::regression_tests::SpreadPointsFromRoot(points, rank);
}

void Build(VoronoiGrid &grid, const std::vector<Vector3D> &points)
{
#ifdef MADVORO_WITH_MPI
    grid.BuildParallel(points);
#else
    grid.Build(points);
#endif
}

struct Errors
{
    std::size_t cells = 0;
    std::size_t badCubes = 0; // cells that are not the cube around their point
    std::size_t unmatchedCells = 0; // cells the Delaunay builder does not have, or has with other neighbors
    double maxVolumeError = 0;
    double maxCMError = 0;
    double totalVolume = 0;
};

// the cell's neighbors over faces that are not tiny
std::vector<Key> NeighborsOf(const VoronoiGrid &grid, std::size_t cell, double h)
{
    std::vector<Key> neighbors;
    for(std::size_t face : grid.GetCellFaces(cell))
    {
        const auto pair = grid.GetFaceNeighbors(face);
        const std::size_t other = (pair.first == cell) ? pair.second : pair.first;
        if(grid.GetArea(face) > TINY_FACE * h * h)
        {
            neighbors.push_back(KeyOf(grid.GetMeshPoint(other)));
        }
    }
    std::sort(neighbors.begin(), neighbors.end());
    return neighbors;
}

// whether the cell has exactly the faces of its cube: 6 neighbors, one at distance h along each axis direction
bool IsCube(const VoronoiGrid &grid, std::size_t cell, double h)
{
    const Vector3D p = grid.GetMeshPoint(cell);
    const std::vector<Key> neighbors = NeighborsOf(grid, cell, h);
    if(neighbors.size() != 6)
    {
        return false;
    }
    int directions = 0;
    for(const Key &n : neighbors)
    {
        const double d[3] = {n[0] - p.x, n[1] - p.y, n[2] - p.z};
        for(int axis = 0; axis < 3; ++axis)
        {
            const double along = d[axis], across = std::abs(d[(axis + 1) % 3]) + std::abs(d[(axis + 2) % 3]);
            if(std::abs(std::abs(along) - h) < 1e-9 * h and across < 1e-9 * h)
            {
                directions |= 1 << (2 * axis + (along > 0 ? 1 : 0));
            }
        }
    }
    return directions == 63;
}

Errors Check(const VoronoiGrid &clipped, const VoronoiGrid &delaunay, double h)
{
    Errors errors;
    std::map<Key, std::size_t> cellOf;
    for(std::size_t c = 0; c < delaunay.GetPointNo(); ++c)
    {
        cellOf[KeyOf(delaunay.GetMeshPoint(c))] = c;
    }
    const double volume = h * h * h;
    errors.cells = clipped.GetPointNo();
    errors.unmatchedCells = (clipped.GetPointNo() != delaunay.GetPointNo()) ? 1 : 0;
    for(std::size_t c = 0; c < clipped.GetPointNo(); ++c)
    {
        errors.totalVolume += clipped.GetVolume(c);
        errors.maxVolumeError = std::max(errors.maxVolumeError, std::abs(clipped.GetVolume(c) - volume) / volume);
        const Vector3D offset = clipped.GetCellCM(c) - clipped.GetMeshPoint(c);
        errors.maxCMError = std::max(errors.maxCMError, std::sqrt(offset.x * offset.x + offset.y * offset.y + offset.z * offset.z) / h);
        errors.badCubes += IsCube(clipped, c, h) ? 0 : 1;
        const auto it = cellOf.find(KeyOf(clipped.GetMeshPoint(c)));
        if(it == cellOf.end() or NeighborsOf(clipped, c, h) != NeighborsOf(delaunay, it->second, h))
        {
            ++errors.unmatchedCells;
        }
    }
    return errors;
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
    // a spacing that is not a power of 2 (no binary fraction is exact)
    const std::size_t side = (worldSize > 1) ? 40 : 20;
    const double h = 1.0 / static_cast<double>(side);
    std::ofstream metrics;
    if(rank == 0)
    {
        metrics.open("voronoi_clipping_cartesian_metrics.txt");
        metrics << "mode " << ((worldSize > 1) ? "mpi" : "serial") << "\n";
    }
    bool pass = true;
    const std::vector<Vector3D> points = MakeGrid(side, rank);
    for(const std::string setup : {"box", "periodic"})
    {
        VoronoiGrid delaunay(Vector3D(0, 0, 0), Vector3D(1, 1, 1)), clipped(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
        clipped.SetCellBuilder(VoronoiGrid::CellBuilder::Clipping);
        for(VoronoiGrid *grid : {&delaunay, &clipped})
        {
            if(setup == "periodic")
            {
                grid->SetPeriodicBoundaries(true, true, true);
            }
        }
        // a first build, then rebuilds of the same points (with the radii and ghosts of the build before)
        for(int build = 0; build < BUILDS; ++build)
        {
            for(VoronoiGrid *grid : {&delaunay, &clipped})
            {
                Build(*grid, (build == 0) ? points : grid->getAllPoints());
            }
            const Errors local = Check(clipped, delaunay, h);
            const std::size_t badCubes = GlobalSum(local.badCubes);
            const std::size_t unmatched = GlobalSum(local.unmatchedCells);
            const double volumeError = GlobalMax(local.maxVolumeError);
            const double cmError = GlobalMax(local.maxCMError);
            double total = local.totalVolume;
#ifdef MADVORO_WITH_MPI
            MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
            const bool ok = badCubes == 0 and unmatched == 0 and volumeError < MOMENTS_TOLERANCE and cmError < MOMENTS_TOLERANCE and std::abs(total - 1.0) < 1e-12;
            pass = pass and ok;
            const std::string name = setup + "_build" + std::to_string(build);
            if(rank == 0)
            {
                std::cout << name << ": bad_cubes=" << badCubes << " unmatched_cells=" << unmatched << " max_volume_error=" << volumeError << " max_cm_error=" << cmError << " total_volume_error=" << std::abs(total - 1.0) << (ok ? "" : " FAIL") << std::endl;
                metrics << name << "_pass " << (ok ? 1 : 0) << "\n";
                if(not ok)
                {
                    metrics << "failed " << name << "\n";
                }
            }
        }
    }

    const int passed = pass ? 1 : 0;
    if(rank == 0)
    {
        metrics << "pass " << passed << "\n";
        std::cout << "voronoi_clipping_cartesian " << (passed ? "PASS" : "FAIL") << std::endl;
    }
#ifdef MADVORO_WITH_MPI
    MPI_Finalize();
#endif
    return passed ? 0 : 1;
}
