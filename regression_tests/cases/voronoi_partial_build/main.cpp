// Partial builds must match a full build on the participating cells and keep the other points.
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "regression_tests/lib/voronoi_test_common.hpp"

#ifdef MADVORO_WITH_MPI
#include <mpi.h>
#endif

namespace
{

using MadVoro::regression_tests::VoronoiGrid;
using PointKey = std::tuple<double, double, double>;

// Gathers (x, y, z, volume) of every owned cell to rank 0.
std::vector<double> GatherCells(const VoronoiGrid &grid)
{
    std::vector<double> local;
    const std::size_t cells = grid.GetPointNo();
    local.reserve(4 * cells);
    for(std::size_t i = 0; i < cells; ++i)
    {
        const Vector3D &point = grid.GetMeshPoint(i);
        local.insert(local.end(), {point.x, point.y, point.z, grid.GetVolume(i)});
    }
#ifdef MADVORO_WITH_MPI
    int rank = 0;
    int worldSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    int localCount = static_cast<int>(local.size());
    std::vector<int> counts(worldSize), displs(worldSize, 0);
    MPI_Gather(&localCount, 1, MPI_INT, counts.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<double> global;
    if(rank == 0)
    {
        for(int r = 1; r < worldSize; ++r)
        {
            displs[r] = displs[r - 1] + counts[r - 1];
        }
        global.resize(displs[worldSize - 1] + counts[worldSize - 1]);
    }
    MPI_Gatherv(local.data(), localCount, MPI_DOUBLE, global.data(), counts.data(), displs.data(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    return global;
#else
    return local;
#endif
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

void Build(VoronoiGrid &grid, const std::vector<Vector3D> &points, const std::vector<std::size_t> &active)
{
#ifdef MADVORO_WITH_MPI
    grid.BuildPartiallyParallel(points, active);
#else
    grid.BuildPartially(points, active);
#endif
}

// Rank- and order-independent subset rules, so every rank selects consistently.
bool InHalfBox(const Vector3D &point)
{
    return point.x < 0.5;
}

bool InScatteredThird(const Vector3D &point)
{
    const std::uint64_t qx = static_cast<std::uint64_t>(point.x * 1048576.0);
    const std::uint64_t qy = static_cast<std::uint64_t>(point.y * 1048576.0);
    const std::uint64_t qz = static_cast<std::uint64_t>(point.z * 1048576.0);
    return ((qx * 73856093ULL) ^ (qy * 19349663ULL) ^ (qz * 83492791ULL)) % 3 == 0;
}

struct SubsetResult
{
    std::size_t activeCells = 0;
    std::size_t expectedActive = 0;
    std::size_t allPointsAfter = 0;
    std::size_t unmatched = 0;
    std::size_t duplicates = 0;
    double maxRelError = 0.0;
};

// Compares one partial build with the full-build volumes; meaningful on rank 0 only.
SubsetResult CheckSubset(VoronoiGrid &grid, bool (*isActive)(const Vector3D &), const std::map<PointKey, double> &fullVolumes)
{
    const std::vector<Vector3D> points = grid.getAllPoints();
    std::vector<std::size_t> active;
    for(std::size_t i = 0; i < points.size(); ++i)
    {
        if(isActive(points[i]))
        {
            active.push_back(i);
        }
    }

    SubsetResult result;
    result.expectedActive = GlobalSum(active.size());
    Build(grid, points, active);
    result.allPointsAfter = GlobalSum(grid.getAllPoints().size());

    const std::vector<double> cells = GatherCells(grid);
    result.activeCells = cells.size() / 4;
    std::map<PointKey, int> seen;
    for(std::size_t c = 0; c + 3 < cells.size(); c += 4)
    {
        const PointKey key(cells[c], cells[c + 1], cells[c + 2]);
        if(++seen[key] > 1)
        {
            result.duplicates++;
        }
        const auto it = fullVolumes.find(key);
        if(it == fullVolumes.end())
        {
            result.unmatched++;
            continue;
        }
        result.maxRelError = std::max(result.maxRelError, std::abs(cells[c + 3] - it->second) / it->second);
    }
    return result;
}

bool Report(int rank, const std::string &label, const SubsetResult &result, std::size_t totalPoints, double tolerance, std::ostream &metrics)
{
    const bool pass = result.activeCells == result.expectedActive and result.expectedActive > 0 and result.expectedActive < totalPoints and result.allPointsAfter == totalPoints and result.unmatched == 0 and result.duplicates == 0 and result.maxRelError < tolerance;
    if(rank == 0)
    {
        std::cout << label << ": active_cells=" << result.activeCells << " expected=" << result.expectedActive << " all_points_after=" << result.allPointsAfter << " (of " << totalPoints << ")" << " unmatched=" << result.unmatched << " duplicates=" << result.duplicates << " max_rel_error=" << result.maxRelError << " pass=" << (pass ? 1 : 0) << std::endl;
        metrics << label << "_max_rel_error " << result.maxRelError << "\n" << label << "_pass " << (pass ? 1 : 0) << "\n";
    }
    return pass;
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
    const std::uint64_t seed = 13579ULL;
    const std::size_t numPoints = (worldSize > 1) ? static_cast<std::size_t>(50000) : static_cast<std::size_t>(20000);
    const double tolerance = 1e-10;

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

    // Full build: the reference volumes. Later builds reuse this grid, as RICH does.
    VoronoiGrid grid(ll, ur);
    MadVoro::regression_tests::BuildVoronoiMesh(grid, points);
    std::map<PointKey, double> fullVolumes;
    {
        const std::vector<double> cells = GatherCells(grid);
        for(std::size_t c = 0; c + 3 < cells.size(); c += 4)
        {
            fullVolumes[PointKey(cells[c], cells[c + 1], cells[c + 2])] = cells[c + 3];
        }
    }

    std::ofstream metrics;
    if(rank == 0)
    {
        metrics.open("voronoi_partial_build_metrics.txt");
        metrics << "mode " << (worldSize > 1 ? "mpi" : "serial") << "\n" << "num_points " << numPoints << "\n";
    }

    bool pass = (rank != 0) or (fullVolumes.size() == numPoints);
    pass = Report(rank, "half_box", CheckSubset(grid, InHalfBox, fullVolumes), numPoints, tolerance, metrics) and pass;
    pass = Report(rank, "scattered_third", CheckSubset(grid, InScatteredThird, fullVolumes), numPoints, tolerance, metrics) and pass;

    // A full rebuild after partial builds must still tile the box.
    {
        const std::vector<Vector3D> current = grid.getAllPoints();
        MadVoro::regression_tests::BuildVoronoiMesh(grid, current);
        const double totalVolume = MadVoro::regression_tests::ReduceGlobalVolume(MadVoro::regression_tests::SumOwnedCellVolumes(grid));
        const double boxVolume = MadVoro::regression_tests::BoxVolume(ll, ur);
        const double relError = std::abs(totalVolume - boxVolume) / boxVolume;
        const std::size_t cells = GlobalSum(grid.GetPointNo());
        const bool rebuildPass = relError < tolerance and cells == numPoints;
        if(rank == 0)
        {
            std::cout << "full_rebuild: cells=" << cells << " rel_volume_error=" << relError << " pass=" << (rebuildPass ? 1 : 0) << std::endl;
            metrics << "full_rebuild_rel_error " << relError << "\n" << "full_rebuild_pass " << (rebuildPass ? 1 : 0) << "\n";
        }
        pass = rebuildPass and pass;
    }

    int passed = pass ? 1 : 0;
#ifdef MADVORO_WITH_MPI
    MPI_Bcast(&passed, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    if(rank == 0)
    {
        metrics << "pass " << passed << "\n";
        metrics.close();
        std::cout << "voronoi_partial_build " << (passed ? "PASS" : "FAIL") << " seed=" << seed << std::endl;
    }

#ifdef MADVORO_WITH_MPI
    MPI_Finalize();
#endif
    return passed ? 0 : 1;
}
