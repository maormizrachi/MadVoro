// On several ranks, builds from thread subdomains must match builds as one domain: same owned points,
// neighbors, close moments, and correct ghost records and exchanged data. Also checked: bit for bit results
// (one subdomain, repeated builds, mixed requests), load balance switches, periodic boxes and partial builds.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <random>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include <mpi.h>
#include <mpi_utils/mpi_exchange.hpp>

#include "regression_tests/lib/voronoi_test_common.hpp"

namespace
{

using MadVoro::regression_tests::VoronoiGrid;
using Key = std::tuple<double, double, double>;

// Relative to the cell's width squared
constexpr double CELL_FACE_AREA = 1e-8;
constexpr double SLIVER_AREA = 1e-6;
constexpr double MOMENTS_TOLERANCE = 1e-7;

Key KeyOf(const Vector3D &p)
{
    return Key(p.x, p.y, p.z);
}

double Distance(const Vector3D &a, const Vector3D &b)
{
    return std::sqrt((a.x - b.x) * (a.x - b.x) + (a.y - b.y) * (a.y - b.y) + (a.z - b.z) * (a.z - b.z));
}

std::vector<Vector3D> MakePoints(const std::string &kind, std::size_t count, int rank)
{
    std::vector<Vector3D> points;
    if(rank == 0)
    {
        if(kind == "lattice")
        {
            // Dyadic lattice: cospherical degeneracy everywhere
            const std::size_t side = 24;
            const double h = 1.0 / 32.0;
            for(std::size_t i = 0; i < side; ++i)
            {
                for(std::size_t j = 0; j < side; ++j)
                {
                    for(std::size_t k = 0; k < side; ++k)
                    {
                        points.emplace_back(h * (i + 4), h * (j + 4), h * (k + 4));
                    }
                }
            }
        }
        else if(kind == "clustered")
        {
            // Most points in a ball, so the ranks' domains differ much in size
            std::mt19937_64 gen(97531ULL);
            std::normal_distribution<double> blob(0.0, 0.08);
            std::uniform_real_distribution<double> uniform(0.0, 1.0);
            const auto clip = [](double x) { return std::min(std::max(x, 1e-6), 1.0 - 1e-6); };
            for(std::size_t i = 0; i < count; ++i)
            {
                if(i % 10 < 7)
                {
                    points.emplace_back(clip(0.3 + blob(gen)), clip(0.35 + blob(gen)), clip(0.4 + blob(gen)));
                }
                else
                {
                    points.emplace_back(uniform(gen), uniform(gen), uniform(gen));
                }
            }
        }
        else
        {
            points = MadVoro::regression_tests::RandRectangular(count, Vector3D(0, 0, 0), Vector3D(1, 1, 1), 13579ULL);
        }
    }
    return MadVoro::regression_tests::SpreadPointsFromRoot(points, rank);
}

// Small deterministic motion (the lattice stays fixed)
std::vector<Vector3D> Moved(const std::vector<Vector3D> &points, const std::string &kind)
{
    if(kind == "lattice")
    {
        return points;
    }
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

void BuildTwice(VoronoiGrid &grid, const std::vector<Vector3D> &points, const std::string &kind)
{
    const std::vector<Vector3D> mine = grid.BuildParallel(points);
    grid.BuildParallel(Moved(mine, kind));
}

struct Setup
{
    std::string name;
    std::array<bool, 3> periodic;
    bool partial;
};

// A scattered third of the points, by position
std::vector<std::size_t> InScatteredThird(const std::vector<Vector3D> &points)
{
    std::vector<std::size_t> active;
    for(std::size_t i = 0; i < points.size(); ++i)
    {
        const std::uint64_t qx = static_cast<std::uint64_t>(points[i].x * 1048576.0);
        const std::uint64_t qy = static_cast<std::uint64_t>(points[i].y * 1048576.0);
        const std::uint64_t qz = static_cast<std::uint64_t>(points[i].z * 1048576.0);
        if(((qx * 73856093ULL) ^ (qy * 19349663ULL) ^ (qz * 83492791ULL)) % 3 == 0)
        {
            active.push_back(i);
        }
    }
    return active;
}

void BuildTwice(VoronoiGrid &grid, const std::vector<Vector3D> &points, const std::string &kind, const Setup &setup)
{
    grid.SetPeriodicBoundaries(setup.periodic[0], setup.periodic[1], setup.periodic[2]);
    if(not setup.partial)
    {
        BuildTwice(grid, points, kind);
        return;
    }
    grid.BuildPartiallyParallel(points, InScatteredThird(points));
    const std::vector<Vector3D> moved = Moved(grid.getAllPoints(), kind);
    grid.BuildPartiallyParallel(moved, InScatteredThird(moved));
}

// Uneven weights by position: the corner near the origin is 20 times heavier
std::vector<double> UnevenWeights(const std::vector<Vector3D> &points)
{
    std::vector<double> weights;
    for(const Vector3D &p : points)
    {
        weights.push_back((p.x < 0.4 and p.y < 0.5) ? 20.0 : 1.0);
    }
    return weights;
}

// Rebalance to uneven weights, switch back (both rebuild with MockMesh), then build on moved points
constexpr std::size_t SWITCH_STEPS = 3;
const char *const SWITCH_STEP_NAMES[SWITCH_STEPS] = {"rebalanced", "switched_back", "rebuilt"};

struct LoadBalanceSwitches
{
    VoronoiGrid grid{Vector3D(0, 0, 0), Vector3D(1, 1, 1)};
    std::shared_ptr<LoadBalancer<Vector3D>> before;
    bool rebalanced = false;

    void Step(std::size_t step, const std::string &kind)
    {
        if(step == 0)
        {
            this->before = this->grid.GetLoadBalancer()->clone();
            this->grid.Rebalance(UnevenWeights(this->grid.getAllPoints()));
            this->rebalanced = this->grid.DidRebalance();
        }
        else if(step == 1)
        {
            this->grid.SetLoadBalancer(this->before);
        }
        else
        {
            this->grid.BuildParallel(Moved(this->grid.getAllPoints(), kind));
        }
    }
};

// Neighbors over non-negligible faces, by position
std::set<Key> Neighbors(const VoronoiGrid &grid, std::size_t cell, double width)
{
    std::set<Key> neighbors;
    for(std::size_t face : grid.GetCellFaces(cell))
    {
        const auto &pair = grid.GetFaceNeighbors(face);
        const std::size_t other = (pair.first == cell) ? pair.second : pair.first;
        const bool mirror = (other >= grid.GetPointNo()) and grid.IsPointOutsideBox(other);
        const double minArea = (mirror ? SLIVER_AREA : CELL_FACE_AREA) * width * width;
        if(grid.GetArea(face) > minArea)
        {
            neighbors.insert(KeyOf(grid.GetMeshPoint(other)));
        }
    }
    return neighbors;
}

struct Comparison
{
    std::size_t ownershipMismatches = 0;
    std::size_t neighborMismatches = 0;
    std::size_t unlistedFaces = 0;
    std::size_t unrecordedGhosts = 0;
    std::size_t wrongExchangedData = 0;
    double maxVolumeError = 0;
    double maxCMError = 0;
    double maxGhostVolumeError = 0;
    double maxGhostCMError = 0;
};

// Neighboring ghosts must be recorded, and data exchanged through the records must land on the right points
void CheckGhosts(VoronoiGrid &grid, Comparison &result)
{
    const std::size_t cells = grid.GetPointNo();
    std::set<std::size_t> recorded;
    for(const std::vector<std::size_t> &indices : grid.GetGhostIndeces())
    {
        recorded.insert(indices.begin(), indices.end());
    }
    // The rank's own points among the ghosts need no exchange
    const auto &local = grid.GetIndicesInAllPoints();
    for(std::size_t c = 0; c < cells; ++c)
    {
        for(std::size_t face : grid.GetCellFaces(c))
        {
            const auto &pair = grid.GetFaceNeighbors(face);
            const std::size_t other = (pair.first == c) ? pair.second : pair.first;
            if(other >= cells and not grid.IsPointOutsideBox(other) and recorded.count(other) == 0 and local.find(other) == local.end())
            {
                result.unrecordedGhosts++;
            }
        }
    }
    // Each rank sends the positions of all its points (not only its cells in a partial build)
    std::vector<double> coordinates;
    for(const Vector3D &p : grid.getAllPoints())
    {
        coordinates.insert(coordinates.end(), {p.x, p.y, p.z});
    }
    const std::vector<std::vector<double>> received = MPI_exchange_data_indexed(grid.GetDuplicatedProcs(), coordinates, grid.GetDuplicatedPoints(), 3);
    const std::vector<std::vector<std::size_t>> &ghosts = grid.GetGhostIndeces();
    for(std::size_t i = 0; i < ghosts.size(); ++i)
    {
        if(received[i].size() != 3 * ghosts[i].size())
        {
            result.wrongExchangedData += 1 + ghosts[i].size();
            continue;
        }
        for(std::size_t j = 0; j < ghosts[i].size(); ++j)
        {
            // A periodic image is its point moved by a period
            const Vector3D &p = grid.GetMeshPoint(ghosts[i][j]);
            const Vector3D &t = grid.GetPeriodicImageTranslation(ghosts[i][j]);
            if(received[i][3 * j] + t.x != p.x or received[i][3 * j + 1] + t.y != p.y or received[i][3 * j + 2] + t.z != p.z)
            {
                result.wrongExchangedData++;
            }
        }
    }
}

Comparison Compare(VoronoiGrid &subdomains, VoronoiGrid &reference)
{
    Comparison result;
    const std::size_t cells = reference.GetPointNo();
    if(subdomains.GetPointNo() != cells)
    {
        result.ownershipMismatches = cells + 1;
        return result;
    }
    std::map<Key, std::size_t> referenceCell;
    for(std::size_t c = 0; c < cells; ++c)
    {
        referenceCell[KeyOf(reference.GetMeshPoint(c))] = c;
    }
    for(std::size_t c = 0; c < cells; ++c)
    {
        const auto it = referenceCell.find(KeyOf(subdomains.GetMeshPoint(c)));
        if(it == referenceCell.end())
        {
            result.ownershipMismatches++;
            continue;
        }
        const std::size_t r = it->second;
        const double width = std::cbrt(reference.GetVolume(r));
        result.neighborMismatches += (Neighbors(subdomains, c, width) != Neighbors(reference, r, width)) ? 1 : 0;
        result.maxVolumeError = std::max(result.maxVolumeError, std::abs(subdomains.GetVolume(c) - reference.GetVolume(r)) / reference.GetVolume(r));
        result.maxCMError = std::max(result.maxCMError, Distance(subdomains.GetCellCM(c), reference.GetCellCM(r)) / width);
    }
    for(std::size_t f = 0; f < subdomains.GetTotalFacesNumber(); ++f)
    {
        const auto &pair = subdomains.GetFaceNeighbors(f);
        for(std::size_t cell : {pair.first, pair.second})
        {
            if(cell >= cells)
            {
                continue;
            }
            const auto &faces = subdomains.GetCellFaces(cell);
            result.unlistedFaces += (std::find(faces.begin(), faces.end(), f) == faces.end()) ? 1 : 0;
        }
    }
    // Ghost moments must be those of their cells on the owning ranks
    std::map<Key, std::pair<double, Vector3D>> referenceGhost;
    for(const std::vector<std::size_t> &indices : reference.GetGhostIndeces())
    {
        for(std::size_t g : indices)
        {
            referenceGhost[KeyOf(reference.GetMeshPoint(g))] = {reference.GetVolume(g), reference.GetCellCM(g)};
        }
    }
    for(const std::vector<std::size_t> &indices : subdomains.GetGhostIndeces())
    {
        for(std::size_t g : indices)
        {
            const auto it = referenceGhost.find(KeyOf(subdomains.GetMeshPoint(g)));
            if(it == referenceGhost.end())
            {
                continue;
            }
            // An unbuilt point in a partial build has volume 0
            const auto &[volume, CM] = it->second;
            const double scale = (volume > 0) ? volume : 1.0;
            result.maxGhostVolumeError = std::max(result.maxGhostVolumeError, std::abs(subdomains.GetVolume(g) - volume) / scale);
            result.maxGhostCMError = std::max(result.maxGhostCMError, Distance(subdomains.GetCellCM(g), CM) / std::cbrt(scale));
        }
    }
    CheckGhosts(subdomains, result);
    return result;
}

bool Same(const Vector3D &a, const Vector3D &b)
{
    return a.x == b.x and a.y == b.y and a.z == b.z;
}

// Differences between two builds, by category
struct Differences
{
    std::map<std::string, std::size_t> counts;

    void Add(const char *what, bool different)
    {
        this->counts[what] += different ? 1 : 0;
    }

    std::size_t Total() const
    {
        std::size_t total = 0;
        for(const auto &entry : this->counts)
        {
            total += entry.second;
        }
        return total;
    }
};

// Differences between two builds that must be the same bit for bit
Differences CountExactDifferences(VoronoiGrid &a, VoronoiGrid &b)
{
    Differences differences;
    const auto count = [&differences](bool different, const char *what) { differences.Add(what, different); };
    count(a.GetPointNo() != b.GetPointNo(), "cells");
    count(a.getMeshPoints().size() != b.getMeshPoints().size(), "mesh_points");
    count(a.GetTotalFacesNumber() != b.GetTotalFacesNumber(), "faces");
    count(a.GetFacePoints().size() != b.GetFacePoints().size(), "face_vertices");
    count(a.GetGhostIndeces() != b.GetGhostIndeces(), "ghost_records");
    count(a.GetDuplicatedProcs() != b.GetDuplicatedProcs(), "sent_procs");
    count(a.GetDuplicatedPoints() != b.GetDuplicatedPoints(), "sent_points");
    if(differences.Total() > 0)
    {
        return differences;
    }
    for(std::size_t p = 0; p < a.getMeshPoints().size(); ++p)
    {
        count(not Same(a.getMeshPoints()[p], b.getMeshPoints()[p]), "mesh_points");
    }
    for(std::size_t c = 0; c < a.GetPointNo(); ++c)
    {
        count(a.GetVolume(c) != b.GetVolume(c), "volumes");
        count(not Same(a.GetCellCM(c), b.GetCellCM(c)), "centroids");
        const auto &facesA = a.GetCellFaces(c);
        const auto &facesB = b.GetCellFaces(c);
        count(not std::equal(facesA.begin(), facesA.end(), facesB.begin(), facesB.end()), "cell_faces");
    }
    for(std::size_t f = 0; f < a.GetTotalFacesNumber(); ++f)
    {
        count(a.GetFaceNeighbors(f) != b.GetFaceNeighbors(f), "face_neighbors");
        count(a.GetArea(f) != b.GetArea(f), "face_areas");
        count(not Same(a.FaceCM(f), b.FaceCM(f)), "face_centroids");
        const auto &verticesA = a.GetPointsInFace(f);
        const auto &verticesB = b.GetPointsInFace(f);
        count(verticesA.size() != verticesB.size(), "face_vertices");
        for(std::size_t k = 0; k < std::min(verticesA.size(), verticesB.size()); ++k)
        {
            count(not Same(a.GetFacePoints()[verticesA[k]], b.GetFacePoints()[verticesB[k]]), "face_vertices");
        }
    }
    return differences;
}

// Differences between two builds that must be the same bit for bit up to the numbering of the ghosts
Differences CountDifferencesUpToGhostNumbering(VoronoiGrid &a, VoronoiGrid &b)
{
    Differences differences;
    const auto count = [&differences](bool different, const char *what) { differences.Add(what, different); };
    count(a.GetPointNo() != b.GetPointNo(), "cells");
    count(a.getMeshPoints().size() != b.getMeshPoints().size(), "mesh_points");
    count(a.GetTotalFacesNumber() != b.GetTotalFacesNumber(), "faces");
    count(a.GetDuplicatedProcs() != b.GetDuplicatedProcs(), "sent_procs");
    count(a.GetDuplicatedPoints() != b.GetDuplicatedPoints(), "sent_points");
    count(a.GetGhostIndeces().size() != b.GetGhostIndeces().size(), "ghost_records");
    if(differences.Total() > 0)
    {
        return differences;
    }
    using Face = std::tuple<Key, double, Key>;
    const auto faces = [](VoronoiGrid &grid, std::size_t cell)
    {
        std::multiset<Face> result;
        for(std::size_t face : grid.GetCellFaces(cell))
        {
            const auto &pair = grid.GetFaceNeighbors(face);
            const std::size_t other = (pair.first == cell) ? pair.second : pair.first;
            result.emplace(KeyOf(grid.GetMeshPoint(other)), grid.GetArea(face), KeyOf(grid.FaceCM(face)));
        }
        return result;
    };
    for(std::size_t c = 0; c < a.GetPointNo(); ++c)
    {
        count(not Same(a.GetMeshPoint(c), b.GetMeshPoint(c)), "own_points");
        count(a.GetVolume(c) != b.GetVolume(c), "volumes");
        count(not Same(a.GetCellCM(c), b.GetCellCM(c)), "centroids");
        count(faces(a, c) != faces(b, c), "cell_faces");
    }
    const auto &ghostsA = a.GetGhostIndeces();
    const auto &ghostsB = b.GetGhostIndeces();
    for(std::size_t i = 0; i < ghostsA.size(); ++i)
    {
        count(ghostsA[i].size() != ghostsB[i].size(), "ghost_records");
        for(std::size_t j = 0; j < std::min(ghostsA[i].size(), ghostsB[i].size()); ++j)
        {
            count(not Same(a.GetMeshPoint(ghostsA[i][j]), b.GetMeshPoint(ghostsB[i][j])), "ghost_records");
        }
    }
    return differences;
}

std::size_t GlobalSum(std::size_t value)
{
    unsigned long long sum = value;
    MPI_Allreduce(MPI_IN_PLACE, &sum, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    return static_cast<std::size_t>(sum);
}

double GlobalMax(double value)
{
    MPI_Allreduce(MPI_IN_PLACE, &value, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    return value;
}

} // namespace

int main()
{
    MPI_Init(nullptr, nullptr);
    int rank = 0, worldSize = 1;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    const std::size_t numPoints = 10000 * static_cast<std::size_t>(worldSize);

    std::ofstream metrics;
    if(rank == 0)
    {
        metrics.open("voronoi_thread_subdomains_mpi_metrics.txt");
        metrics << "mode mpi\nranks " << worldSize << "\n";
    }
    bool pass = true;
    for(const std::string kind : {"uniform", "lattice", "clustered"})
    {
        const std::vector<Vector3D> points = MakePoints(kind, numPoints, rank);
        VoronoiGrid reference(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
        BuildTwice(reference, points, kind);
        Comparison referenceGhosts;
        CheckGhosts(reference, referenceGhosts);
        if(GlobalSum(referenceGhosts.unrecordedGhosts + referenceGhosts.wrongExchangedData) != 0 and rank == 0)
        {
            std::cout << kind << ": the build as one domain fails the ghost checks too" << std::endl;
        }

        // Comparison with a build as one domain, within tolerances
        const auto compareWith = [&](VoronoiGrid &grid, VoronoiGrid &oneDomain, const std::string &name)
        {
            const Comparison local = Compare(grid, oneDomain);
            Comparison c;
            c.ownershipMismatches = GlobalSum(local.ownershipMismatches);
            c.neighborMismatches = GlobalSum(local.neighborMismatches);
            c.unlistedFaces = GlobalSum(local.unlistedFaces);
            c.unrecordedGhosts = GlobalSum(local.unrecordedGhosts);
            c.wrongExchangedData = GlobalSum(local.wrongExchangedData);
            c.maxVolumeError = GlobalMax(local.maxVolumeError);
            c.maxCMError = GlobalMax(local.maxCMError);
            c.maxGhostVolumeError = GlobalMax(local.maxGhostVolumeError);
            c.maxGhostCMError = GlobalMax(local.maxGhostCMError);
            double volume = 0, oneDomainVolume = 0;
            for(std::size_t cell = 0; cell < grid.GetPointNo(); ++cell)
            {
                volume += grid.GetVolume(cell);
            }
            for(std::size_t cell = 0; cell < oneDomain.GetPointNo(); ++cell)
            {
                oneDomainVolume += oneDomain.GetVolume(cell);
            }
            oneDomainVolume = MadVoro::regression_tests::ReduceGlobalVolume(oneDomainVolume);
            const double totalVolumeError = std::abs(MadVoro::regression_tests::ReduceGlobalVolume(volume) - oneDomainVolume) / oneDomainVolume;
            const bool ok = c.ownershipMismatches == 0 and c.neighborMismatches == 0 and c.unlistedFaces == 0 and c.unrecordedGhosts == 0 and c.wrongExchangedData == 0 and c.maxVolumeError < MOMENTS_TOLERANCE and c.maxCMError < MOMENTS_TOLERANCE and c.maxGhostVolumeError < MOMENTS_TOLERANCE and c.maxGhostCMError < MOMENTS_TOLERANCE and totalVolumeError < MOMENTS_TOLERANCE;
            pass = pass and ok;
            if(rank == 0)
            {
                std::cout << kind << " " << name << ": ownership_mismatches=" << c.ownershipMismatches << " neighbor_mismatches=" << c.neighborMismatches << " unlisted_faces=" << c.unlistedFaces << " unrecorded_ghosts=" << c.unrecordedGhosts << " wrong_exchanged_data=" << c.wrongExchangedData << " max_volume_error=" << c.maxVolumeError << " max_cm_error=" << c.maxCMError << " max_ghost_volume_error=" << c.maxGhostVolumeError << " max_ghost_cm_error=" << c.maxGhostCMError << " total_volume_error=" << totalVolumeError << (ok ? "" : " FAIL") << std::endl;
                const std::string prefix = kind + "_" + name + "_";
                metrics << prefix << "ownership_mismatches " << c.ownershipMismatches << "\n" << prefix << "neighbor_mismatches " << c.neighborMismatches << "\n" << prefix << "unlisted_faces " << c.unlistedFaces << "\n" << prefix << "unrecorded_ghosts " << c.unrecordedGhosts << "\n" << prefix << "wrong_exchanged_data " << c.wrongExchangedData << "\n" << prefix << "max_volume_error " << c.maxVolumeError << "\n" << prefix << "max_cm_error " << c.maxCMError << "\n" << prefix << "max_ghost_volume_error " << c.maxGhostVolumeError << "\n" << prefix << "max_ghost_cm_error " << c.maxGhostCMError << "\n" << prefix << "total_volume_error " << totalVolumeError << "\n";
                if(not ok)
                {
                    metrics << "failed " << prefix << "comparison\n";
                }
            }
        };
        const auto compare = [&](VoronoiGrid &grid, const std::string &name) { compareWith(grid, reference, name); };
        // Comparison that must hold bit for bit (optionally up to the ghost numbering)
        const auto compareExactly = [&](VoronoiGrid &a, VoronoiGrid &b, const std::string &name, bool sameGhostNumbering)
        {
            const Differences local = sameGhostNumbering ? CountExactDifferences(a, b) : CountDifferencesUpToGhostNumbering(a, b);
            // Every rank reduces the same categories
            std::string what;
            std::size_t differences = 0;
            for(const char *category : {"cells", "mesh_points", "own_points", "faces", "face_vertices", "ghost_records", "sent_procs", "sent_points", "volumes", "centroids", "cell_faces", "face_neighbors", "face_areas", "face_centroids"})
            {
                const auto it = local.counts.find(category);
                const std::size_t global = GlobalSum(it == local.counts.end() ? 0 : it->second);
                differences += global;
                if(global > 0)
                {
                    what += std::string(" ") + category + "=" + std::to_string(global);
                }
            }
            pass = pass and (differences == 0);
            if(rank == 0)
            {
                std::cout << kind << " " << name << ": differences=" << differences << what << (differences == 0 ? "" : " FAIL") << std::endl;
                metrics << kind << "_" << name << "_differences " << differences << "\n";
                if(differences != 0)
                {
                    metrics << "failed " << kind << "_" << name << "\n";
                }
            }
        };

        for(std::size_t count : {1, 2, 4, 7})
        {
            VoronoiGrid grid(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
            grid.SetThreadSubdomains(count);
            BuildTwice(grid, points, kind);
            compare(grid, "subdomains=" + std::to_string(count));
            if(count == 1)
            {
                compareExactly(grid, reference, "subdomains=1_vs_one_domain", false);
            }
            if(count == 4)
            {
                VoronoiGrid again(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
                again.SetThreadSubdomains(count);
                BuildTwice(again, points, kind);
                compareExactly(again, grid, "subdomains=4_repeated", true);
            }
        }

        // Each rank its own number of subdomains
        {
            VoronoiGrid grid(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
            grid.SetThreadSubdomains(1 + static_cast<std::size_t>(rank % 3));
            BuildTwice(grid, points, kind);
            compare(grid, "subdomains=1+rank%3");
        }

        // Some ranks ask for subdomains, the others do not
        if(worldSize > 1)
        {
            VoronoiGrid grid(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
            grid.SetThreadSubdomains((rank % 2 == 1) ? 4 : 0);
            BuildTwice(grid, points, kind);
            compare(grid, "mixed_requests");
        }

        // Load balance switches, each step compared with the same steps on a mesh built as one domain
        {
            std::vector<std::unique_ptr<LoadBalanceSwitches>> switches;
            std::vector<std::string> names;
            const auto add = [&](std::size_t subdomains, const std::string &name)
            {
                switches.push_back(std::make_unique<LoadBalanceSwitches>());
                switches.back()->grid.SetThreadSubdomains(subdomains);
                BuildTwice(switches.back()->grid, points, kind);
                names.push_back(name);
            };
            add(0, "one_domain");
            for(std::size_t count : {1, 2, 4, 7})
            {
                add(count, "subdomains=" + std::to_string(count));
            }
            add(4, "subdomains=4_repeated");
            if(worldSize > 1)
            {
                add((rank % 2 == 1) ? 4 : 0, "mixed_requests");
            }
            for(std::size_t step = 0; step < SWITCH_STEPS; ++step)
            {
                for(std::unique_ptr<LoadBalanceSwitches> &s : switches)
                {
                    s->Step(step, kind);
                }
                VoronoiGrid &oneDomain = switches.front()->grid;
                const std::string at = std::string("switch_") + SWITCH_STEP_NAMES[step] + "_";
                if(step == 0)
                {
                    // The weights must move points between ranks
                    const std::size_t notRebalanced = GlobalSum(switches.front()->rebalanced ? 0 : 1);
                    pass = pass and (notRebalanced == 0);
                    if(rank == 0 and notRebalanced != 0)
                    {
                        std::cout << kind << " " << at << "did_rebalance=0 FAIL" << std::endl;
                        metrics << "failed " << kind << "_" << at << "did_rebalance\n";
                    }
                }
                for(std::size_t i = 1; i < switches.size(); ++i)
                {
                    VoronoiGrid &grid = switches[i]->grid;
                    if(names[i] == "subdomains=4_repeated")
                    {
                        const std::size_t first = std::find(names.begin(), names.end(), "subdomains=4") - names.begin();
                        compareExactly(grid, switches[first]->grid, at + names[i], true);
                    }
                    else
                    {
                        compareWith(grid, oneDomain, at + names[i]);
                    }
                }
            }
        }

        // Periodic boxes and partial builds
        for(const Setup &setup : {Setup{"periodic", {true, true, true}, false}, Setup{"periodic_x", {true, false, false}, false}, Setup{"partial", {false, false, false}, true}, Setup{"partial_periodic", {true, true, true}, true}})
        {
            VoronoiGrid oneDomain(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
            BuildTwice(oneDomain, points, kind, setup);
            Comparison oneDomainGhosts;
            CheckGhosts(oneDomain, oneDomainGhosts);
            if(GlobalSum(oneDomainGhosts.unrecordedGhosts + oneDomainGhosts.wrongExchangedData) != 0 and rank == 0)
            {
                std::cout << kind << " " << setup.name << ": the build as one domain fails the ghost checks too" << std::endl;
            }
            for(std::size_t count : {1, 4, 7})
            {
                VoronoiGrid grid(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
                grid.SetThreadSubdomains(count);
                BuildTwice(grid, points, kind, setup);
                compareWith(grid, oneDomain, setup.name + "_subdomains=" + std::to_string(count));
                if(count == 4)
                {
                    VoronoiGrid again(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
                    again.SetThreadSubdomains(count);
                    BuildTwice(again, points, kind, setup);
                    compareExactly(again, grid, setup.name + "_subdomains=4_repeated", true);
                }
            }
            if(worldSize > 1)
            {
                VoronoiGrid grid(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
                grid.SetThreadSubdomains((rank % 2 == 1) ? 4 : 0);
                BuildTwice(grid, points, kind, setup);
                compareWith(grid, oneDomain, setup.name + "_mixed_requests");
            }
        }
    }

    const int passed = pass ? 1 : 0;
    if(rank == 0)
    {
        metrics << "pass " << passed << "\n";
        std::cout << "voronoi_thread_subdomains_mpi " << (passed ? "PASS" : "FAIL") << std::endl;
    }
    MPI_Finalize();
    return passed ? 0 : 1;
}
