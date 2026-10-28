// Building cells by clipping gives the cells of the Delaunay builder: the same neighbors (over faces that are not tiny), and volumes and centroids to round-off, for several inputs, boxes, builds, thread subdomains and load balance switches
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

#include "regression_tests/lib/voronoi_test_common.hpp"

#ifdef MADVORO_WITH_MPI
#include <mpi.h>
#endif

namespace
{

using MadVoro::regression_tests::VoronoiGrid;
using Key = std::array<double, 3>;

// as between builds of one domain and of thread subdomains (points 1e-9 from a wall make slivers)
constexpr double MOMENTS_TOLERANCE = 1e-7;
// faces this small (relative to their cell's size) may be dropped by one builder only; larger slivers with mirrors
constexpr double TINY_FACE = 1e-8;
constexpr double MIRROR_SLIVER = 1e-6;

Key KeyOf(const Vector3D &p)
{
    return {p.x, p.y, p.z};
}

double Distance(const Vector3D &a, const Vector3D &b)
{
    const Vector3D d = a - b;
    return std::sqrt(d.x * d.x + d.y * d.y + d.z * d.z);
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

std::vector<Vector3D> MakePoints(const std::string &kind, std::size_t count, int rank)
{
    std::vector<Vector3D> points;
    if(rank == 0)
    {
        if(kind == "lattice")
        {
            // a dyadic lattice: cospherical points everywhere
            const std::size_t side = static_cast<std::size_t>(std::cbrt(static_cast<double>(count)));
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
        else if(kind == "clustered")
        {
            // dense clusters in a sparse background: cells of very different sizes side by side
            boost::mt19937_64 gen(24680ULL);
            boost::random::normal_distribution<double> normal(0.0, 0.03);
            points = MadVoro::regression_tests::RandRectangular(count / 4, Vector3D(0, 0, 0), Vector3D(1, 1, 1), 13579ULL);
            const std::vector<Vector3D> centers = {Vector3D(0.3, 0.3, 0.3), Vector3D(0.7, 0.6, 0.4), Vector3D(0.45, 0.8, 0.75)};
            while(points.size() < count)
            {
                const Vector3D &center = centers[points.size() % centers.size()];
                const Vector3D p(center.x + normal(gen), center.y + normal(gen), center.z + normal(gen));
                if(p.x > 0 and p.x < 1 and p.y > 0 and p.y < 1 and p.z > 0 and p.z < 1)
                {
                    points.push_back(p);
                }
            }
        }
        else
        {
            points = MadVoro::regression_tests::RandRectangular(count, Vector3D(0, 0, 0), Vector3D(1, 1, 1), 97531ULL);
        }
    }
#ifdef MADVORO_WITH_MPI
    points = MadVoro::regression_tests::SpreadPointsFromRoot(points, rank);
#endif
    return points;
}

std::vector<Vector3D> Moved(const std::vector<Vector3D> &points)
{
    std::vector<Vector3D> moved = points;
    for(Vector3D &p : moved)
    {
        const double phase = 1e4 * (p.x + 2 * p.y + 3 * p.z);
        p.x = std::min(std::max(p.x + 1e-3 * std::sin(phase), 1e-9), 1.0 - 1e-9);
        p.y = std::min(std::max(p.y + 1e-3 * std::cos(1.3 * phase), 1e-9), 1.0 - 1e-9);
        p.z = std::min(std::max(p.z + 1e-3 * std::sin(0.7 * phase + 1.0), 1e-9), 1.0 - 1e-9);
    }
    return moved;
}

// a rank- and order-independent third of the points
std::vector<std::size_t> Third(const std::vector<Vector3D> &points)
{
    std::vector<std::size_t> active;
    for(std::size_t i = 0; i < points.size(); ++i)
    {
        const unsigned long long qx = static_cast<unsigned long long>(points[i].x * 1048576.0);
        const unsigned long long qy = static_cast<unsigned long long>(points[i].y * 1048576.0);
        const unsigned long long qz = static_cast<unsigned long long>(points[i].z * 1048576.0);
        if(((qx * 73856093ULL) ^ (qy * 19349663ULL) ^ (qz * 83492791ULL)) % 3 == 0)
        {
            active.push_back(i);
        }
    }
    return active;
}

void Build(VoronoiGrid &grid, const std::vector<Vector3D> &points, bool partial)
{
#ifdef MADVORO_WITH_MPI
    if(partial)
    {
        grid.BuildPartiallyParallel(points, Third(points));
    }
    else
    {
        grid.BuildParallel(points);
    }
#else
    if(partial)
    {
        grid.BuildPartially(points, Third(points));
    }
    else
    {
        grid.Build(points);
    }
#endif
}

std::set<Key> NeighborsOf(const VoronoiGrid &grid, std::size_t cell)
{
    std::set<Key> neighbors;
    const double scale = std::pow(grid.GetVolume(cell), 2.0 / 3.0);
    for(std::size_t face : grid.GetCellFaces(cell))
    {
        const auto pair = grid.GetFaceNeighbors(face);
        const std::size_t other = (pair.first == cell) ? pair.second : pair.first;
        const bool mirror = (other >= grid.GetPointNo()) and grid.IsPointOutsideBox(other);
        if(grid.GetArea(face) > (mirror ? MIRROR_SLIVER : TINY_FACE) * scale)
        {
            neighbors.insert(KeyOf(grid.GetMeshPoint(other)));
        }
    }
    return neighbors;
}

struct Comparison
{
    std::size_t cells = 0;
    std::size_t unmatchedCells = 0;
    std::size_t neighborMismatches = 0;
    double maxVolumeError = 0;
    double maxCMError = 0;
    double totalVolume = 0;
};

Comparison Compare(const VoronoiGrid &clipped, const VoronoiGrid &delaunay)
{
    Comparison result;
    std::map<Key, std::size_t> cellOf;
    for(std::size_t c = 0; c < delaunay.GetPointNo(); ++c)
    {
        cellOf[KeyOf(delaunay.GetMeshPoint(c))] = c;
    }
    result.cells = clipped.GetPointNo();
    result.unmatchedCells = (clipped.GetPointNo() != delaunay.GetPointNo()) ? 1 : 0;
    for(std::size_t c = 0; c < clipped.GetPointNo(); ++c)
    {
        result.totalVolume += clipped.GetVolume(c);
        const auto it = cellOf.find(KeyOf(clipped.GetMeshPoint(c)));
        if(it == cellOf.end())
        {
            ++result.unmatchedCells;
            continue;
        }
        const std::size_t d = it->second;
        const double volume = delaunay.GetVolume(d);
        result.maxVolumeError = std::max(result.maxVolumeError, std::abs(clipped.GetVolume(c) - volume) / volume);
        result.maxCMError = std::max(result.maxCMError, Distance(clipped.GetCellCM(c), delaunay.GetCellCM(d)) / std::cbrt(volume));
        // a neighbor across a face that is tiny in one of the two may be missing from the other
        const std::set<Key> a = NeighborsOf(clipped, c), b = NeighborsOf(delaunay, d);
        for(const Key &k : a)
        {
            if(b.count(k) == 0)
            {
                bool tinyInOther = false;
                for(std::size_t face : delaunay.GetCellFaces(d))
                {
                    const auto pair = delaunay.GetFaceNeighbors(face);
                    tinyInOther = tinyInOther or KeyOf(delaunay.GetMeshPoint((pair.first == d) ? pair.second : pair.first)) == k;
                }
                result.neighborMismatches += tinyInOther ? 0 : 1;
            }
        }
        for(const Key &k : b)
        {
            if(a.count(k) == 0)
            {
                bool tinyInOther = false;
                for(std::size_t face : clipped.GetCellFaces(c))
                {
                    const auto pair = clipped.GetFaceNeighbors(face);
                    tinyInOther = tinyInOther or KeyOf(clipped.GetMeshPoint((pair.first == c) ? pair.second : pair.first)) == k;
                }
                result.neighborMismatches += tinyInOther ? 0 : 1;
            }
        }
    }
    return result;
}

std::vector<double> UnevenWeights(const std::vector<Vector3D> &points)
{
    std::vector<double> weights;
    for(const Vector3D &p : points)
    {
        weights.push_back((p.x < 0.4 and p.y < 0.5) ? 4.0 : 1.0);
    }
    return weights;
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
    const std::size_t numPoints = (worldSize > 1) ? 40000 : 8000;
    std::ofstream metrics;
    if(rank == 0)
    {
        metrics.open("voronoi_clipping_builder_metrics.txt");
        metrics << "mode " << ((worldSize > 1) ? "mpi" : "serial") << "\n";
    }
    bool pass = true;
    const auto report = [&](const std::string &name, const Comparison &local, bool checkTotal)
    {
        const std::size_t unmatched = GlobalSum(local.unmatchedCells);
        const std::size_t mismatches = GlobalSum(local.neighborMismatches);
        const double volumeError = GlobalMax(local.maxVolumeError);
        const double cmError = GlobalMax(local.maxCMError);
        double total = local.totalVolume;
#ifdef MADVORO_WITH_MPI
        MPI_Allreduce(MPI_IN_PLACE, &total, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
#endif
        const double totalError = checkTotal ? std::abs(total - 1.0) : 0.0;
        const bool ok = unmatched == 0 and mismatches == 0 and volumeError < MOMENTS_TOLERANCE and cmError < MOMENTS_TOLERANCE and totalError < 1e-12;
        pass = pass and ok;
        if(rank == 0)
        {
            std::cout << name << ": unmatched_cells=" << unmatched << " neighbor_mismatches=" << mismatches << " max_volume_error=" << volumeError << " max_cm_error=" << cmError << " total_volume_error=" << totalError << (ok ? "" : " FAIL") << std::endl;
            metrics << name << "_pass " << (ok ? 1 : 0) << "\n";
            if(not ok)
            {
                metrics << "failed " << name << "\n";
            }
        }
    };

    for(const std::string kind : {"uniform", "clustered", "lattice"})
    {
        const std::vector<Vector3D> points = MakePoints(kind, numPoints, rank);
        for(const std::string setup : {"box", "periodic", "partial"})
        {
            for(const std::size_t subdomains : {1, 4})
            {
                const std::string name = kind + "_" + setup + "_subdomains=" + std::to_string(subdomains);
                VoronoiGrid delaunay(Vector3D(0, 0, 0), Vector3D(1, 1, 1)), clipped(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
                clipped.SetCellBuilder(VoronoiGrid::CellBuilder::Clipping);
                for(VoronoiGrid *grid : {&delaunay, &clipped})
                {
                    grid->SetThreadSubdomains(subdomains);
                    if(setup == "periodic")
                    {
                        grid->SetPeriodicBoundaries(true, true, true);
                    }
                }
                // twice, the second with the radii and ghosts of the first (the lattice stays degenerate)
                for(VoronoiGrid *grid : {&delaunay, &clipped})
                {
                    Build(*grid, points, setup == "partial");
                }
                for(VoronoiGrid *grid : {&delaunay, &clipped})
                {
                    const std::vector<Vector3D> mine = grid->getAllPoints();
                    Build(*grid, (kind == "lattice") ? mine : Moved(mine), setup == "partial");
                }
                report(name, Compare(clipped, delaunay), setup != "partial");
#ifdef MADVORO_WITH_MPI
                // a load balance switch builds with every ghost known
                if(worldSize > 1 and setup != "partial")
                {
                    for(VoronoiGrid *grid : {&delaunay, &clipped})
                    {
                        grid->Rebalance(UnevenWeights(grid->getAllPoints()));
                    }
                    report(name + "_rebalanced", Compare(clipped, delaunay), true);
                }
#endif
            }
        }
    }

    const int passed = pass ? 1 : 0;
    if(rank == 0)
    {
        metrics << "pass " << passed << "\n";
        std::cout << "voronoi_clipping_builder " << (passed ? "PASS" : "FAIL") << std::endl;
    }
#ifdef MADVORO_WITH_MPI
    MPI_Finalize();
#endif
    return passed ? 0 : 1;
}
