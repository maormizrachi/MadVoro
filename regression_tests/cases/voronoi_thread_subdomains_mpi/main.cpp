// On several ranks, building each rank's mesh from thread subdomains must give the
// tessellation of building each rank as one domain: each rank must own the same points,
// each cell must have the same neighbors and its volume and centroid must agree closely,
// and every face between two cells must be listed by both. The ghosts must be recorded
// for the exchange of cell data between ranks: every point another rank sent that is a
// neighbor of a cell must be among the ghost indices, the data exchanged through the
// ghost indices must be that of the right points, and the ghosts' volumes and centroids
// must be those of their cells on the ranks that own them. Each case builds twice: on the points, then
// on slightly moved points.
// Near a wall a build may leave zero-area sliver faces between a cell and the mirror of
// another point, and these differ between the two builds (see voronoi_thread_subdomains):
// faces against mirrors are compared above SLIVER_AREA.
// Besides that comparison, for uniform, lattice and clustered (uneven between ranks)
// points:
// - one subdomain per rank must give the build as one domain bit for bit, up to the
//   numbering of the ghosts;
// - building again with the same subdomains must give the same result bit for bit (the
//   subdomains run on host threads; a race would show here first);
// - ranks may ask for different numbers of subdomains;
// - when some ranks ask for subdomains and others do not, every rank must build as one
//   domain (the choice is collective), bit for bit.
// The same holds through load balance switches, which rebuild the mesh with MockMesh
// (with subdomains, from them, with the ghosts the old mesh tells it): a rebalance to
// uneven weights, a switch back to the load balance before it, then a build on moved
// points, each compared with the same steps on a mesh built as one domain.
#include <algorithm>
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
            // A dyadic lattice: cospherical degeneracy everywhere
            const std::size_t side = 24;
            const double h = 1.0 / 32.0;
            for(std::size_t i = 0; i < side; ++i)
                for(std::size_t j = 0; j < side; ++j)
                    for(std::size_t k = 0; k < side; ++k)
                        points.emplace_back(h * (i + 4), h * (j + 4), h * (k + 4));
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
                    points.emplace_back(clip(0.3 + blob(gen)), clip(0.35 + blob(gen)), clip(0.4 + blob(gen)));
                else
                    points.emplace_back(uniform(gen), uniform(gen), uniform(gen));
            }
        }
        else
        {
            points = MadVoro::regression_tests::RandRectangular(count, Vector3D(0, 0, 0), Vector3D(1, 1, 1), 13579ULL);
        }
    }
    return MadVoro::regression_tests::SpreadPointsFromRoot(points, rank);
}

// Small deterministic motion of the points a rank holds (the lattice stays fixed)
std::vector<Vector3D> Moved(const std::vector<Vector3D> &points, const std::string &kind)
{
    if(kind == "lattice")
        return points;
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

// Uneven weights by position: the corner near the origin is 20 times heavier
std::vector<double> UnevenWeights(const std::vector<Vector3D> &points)
{
    std::vector<double> weights;
    for(const Vector3D &p : points)
        weights.push_back((p.x < 0.4 and p.y < 0.5) ? 20.0 : 1.0);
    return weights;
}

// The load balance switches, step by step: a rebalance to uneven weights, a switch back to
// the load balance before it (both rebuild the mesh with MockMesh), then a build on moved
// points
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
            this->grid.SetLoadBalancer(this->before);
        else
            this->grid.BuildParallel(Moved(this->grid.getAllPoints(), kind));
    }
};

// A cell's neighbors over faces that are not negligibly small, by position
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
            neighbors.insert(KeyOf(grid.GetMeshPoint(other)));
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

// The ghost bookkeeping of `grid`: the ghosts that are neighbors of cells must be
// recorded, and data exchanged through the records must be that of the right points
void CheckGhosts(VoronoiGrid &grid, Comparison &result)
{
    const std::size_t cells = grid.GetPointNo();
    std::set<std::size_t> recorded;
    for(const std::vector<std::size_t> &indices : grid.GetGhostIndeces())
        recorded.insert(indices.begin(), indices.end());
    for(std::size_t c = 0; c < cells; ++c)
    {
        for(std::size_t face : grid.GetCellFaces(c))
        {
            const auto &pair = grid.GetFaceNeighbors(face);
            const std::size_t other = (pair.first == c) ? pair.second : pair.first;
            if(other >= cells and not grid.IsPointOutsideBox(other) and recorded.count(other) == 0)
                result.unrecordedGhosts++;
        }
    }
    // each rank sends its points' positions; they must land on the ghosts at those positions
    std::vector<double> coordinates;
    for(std::size_t c = 0; c < cells; ++c)
    {
        const Vector3D &p = grid.GetMeshPoint(c);
        coordinates.insert(coordinates.end(), {p.x, p.y, p.z});
    }
    const std::vector<std::vector<double>> received =
        MPI_exchange_data_indexed(grid.GetDuplicatedProcs(), coordinates, grid.GetDuplicatedPoints(), 3);
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
            const Vector3D &p = grid.GetMeshPoint(ghosts[i][j]);
            if(received[i][3 * j] != p.x or received[i][3 * j + 1] != p.y or received[i][3 * j + 2] != p.z)
                result.wrongExchangedData++;
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
        referenceCell[KeyOf(reference.GetMeshPoint(c))] = c;
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
                continue;
            const auto &faces = subdomains.GetCellFaces(cell);
            result.unlistedFaces += (std::find(faces.begin(), faces.end(), f) == faces.end()) ? 1 : 0;
        }
    }
    // ghost volumes and centroids are those of their cells on the ranks that own them
    std::map<Key, std::pair<double, Vector3D>> referenceGhost;
    for(const std::vector<std::size_t> &indices : reference.GetGhostIndeces())
        for(std::size_t g : indices)
            referenceGhost[KeyOf(reference.GetMeshPoint(g))] = {reference.GetVolume(g), reference.GetCellCM(g)};
    for(const std::vector<std::size_t> &indices : subdomains.GetGhostIndeces())
    {
        for(std::size_t g : indices)
        {
            const auto it = referenceGhost.find(KeyOf(subdomains.GetMeshPoint(g)));
            if(it == referenceGhost.end())
                continue;
            const auto &[volume, CM] = it->second;
            result.maxGhostVolumeError = std::max(result.maxGhostVolumeError, std::abs(subdomains.GetVolume(g) - volume) / volume);
            result.maxGhostCMError = std::max(result.maxGhostCMError, Distance(subdomains.GetCellCM(g), CM) / std::cbrt(volume));
        }
    }
    CheckGhosts(subdomains, result);
    return result;
}

bool Same(const Vector3D &a, const Vector3D &b)
{
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

// Differences between two builds, by what differs
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
            total += entry.second;
        return total;
    }
};

// Differences between two builds of a rank that must be the same bit for bit: its
// points, cells, faces, face vertices, moments and the records of its ghosts
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
        return differences;
    for(std::size_t p = 0; p < a.getMeshPoints().size(); ++p)
        count(!Same(a.getMeshPoints()[p], b.getMeshPoints()[p]), "mesh_points");
    for(std::size_t c = 0; c < a.GetPointNo(); ++c)
    {
        count(a.GetVolume(c) != b.GetVolume(c), "volumes");
        count(!Same(a.GetCellCM(c), b.GetCellCM(c)), "centroids");
        const auto &facesA = a.GetCellFaces(c);
        const auto &facesB = b.GetCellFaces(c);
        count(!std::equal(facesA.begin(), facesA.end(), facesB.begin(), facesB.end()), "cell_faces");
    }
    for(std::size_t f = 0; f < a.GetTotalFacesNumber(); ++f)
    {
        count(a.GetFaceNeighbors(f) != b.GetFaceNeighbors(f), "face_neighbors");
        count(a.GetArea(f) != b.GetArea(f), "face_areas");
        count(!Same(a.FaceCM(f), b.FaceCM(f)), "face_centroids");
        const auto &verticesA = a.GetPointsInFace(f);
        const auto &verticesB = b.GetPointsInFace(f);
        count(verticesA.size() != verticesB.size(), "face_vertices");
        for(std::size_t k = 0; k < std::min(verticesA.size(), verticesB.size()); ++k)
            count(!Same(a.GetFacePoints()[verticesA[k]], b.GetFacePoints()[verticesB[k]]), "face_vertices");
    }
    return differences;
}

// Differences between two builds of a rank that must be the same bit for bit up to the
// numbering of the ghosts: a build from subdomains numbers the points other ranks sent
// before the mirrors, a build as one domain in the order it added them. Own points, the
// cells' moments, the faces of each cell (by the position of the point across the face),
// the positions of the ghosts in each record, and the points sent must all be the same.
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
        return differences;
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
        count(!Same(a.GetMeshPoint(c), b.GetMeshPoint(c)), "own_points");
        count(a.GetVolume(c) != b.GetVolume(c), "volumes");
        count(!Same(a.GetCellCM(c), b.GetCellCM(c)), "centroids");
        count(faces(a, c) != faces(b, c), "cell_faces");
    }
    const auto &ghostsA = a.GetGhostIndeces();
    const auto &ghostsB = b.GetGhostIndeces();
    for(std::size_t i = 0; i < ghostsA.size(); ++i)
    {
        count(ghostsA[i].size() != ghostsB[i].size(), "ghost_records");
        for(std::size_t j = 0; j < std::min(ghostsA[i].size(), ghostsB[i].size()); ++j)
            count(!Same(a.GetMeshPoint(ghostsA[i][j]), b.GetMeshPoint(ghostsB[i][j])), "ghost_records");
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
            std::cout << kind << ": the build as one domain fails the ghost checks too" << std::endl;

        // A comparison with a build as one domain, within tolerances
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
            double volume = 0;
            for(std::size_t cell = 0; cell < grid.GetPointNo(); ++cell)
                volume += grid.GetVolume(cell);
            const double totalVolumeError = std::abs(MadVoro::regression_tests::ReduceGlobalVolume(volume) - 1.0);
            const bool ok = c.ownershipMismatches == 0 && c.neighborMismatches == 0 && c.unlistedFaces == 0 && c.unrecordedGhosts == 0 &&
                            c.wrongExchangedData == 0 && c.maxVolumeError < MOMENTS_TOLERANCE && c.maxCMError < MOMENTS_TOLERANCE &&
                            c.maxGhostVolumeError < MOMENTS_TOLERANCE && c.maxGhostCMError < MOMENTS_TOLERANCE &&
                            totalVolumeError < MOMENTS_TOLERANCE;
            pass = pass && ok;
            if(rank == 0)
            {
                std::cout << kind << " " << name << ": ownership_mismatches=" << c.ownershipMismatches
                          << " neighbor_mismatches=" << c.neighborMismatches << " unlisted_faces=" << c.unlistedFaces
                          << " unrecorded_ghosts=" << c.unrecordedGhosts << " wrong_exchanged_data=" << c.wrongExchangedData
                          << " max_volume_error=" << c.maxVolumeError << " max_cm_error=" << c.maxCMError
                          << " max_ghost_volume_error=" << c.maxGhostVolumeError << " max_ghost_cm_error=" << c.maxGhostCMError
                          << " total_volume_error=" << totalVolumeError
                          << (ok ? "" : " FAIL") << std::endl;
                const std::string prefix = kind + "_" + name + "_";
                metrics << prefix << "ownership_mismatches " << c.ownershipMismatches << "\n"
                        << prefix << "neighbor_mismatches " << c.neighborMismatches << "\n"
                        << prefix << "unlisted_faces " << c.unlistedFaces << "\n"
                        << prefix << "unrecorded_ghosts " << c.unrecordedGhosts << "\n"
                        << prefix << "wrong_exchanged_data " << c.wrongExchangedData << "\n"
                        << prefix << "max_volume_error " << c.maxVolumeError << "\n"
                        << prefix << "max_cm_error " << c.maxCMError << "\n"
                        << prefix << "max_ghost_volume_error " << c.maxGhostVolumeError << "\n"
                        << prefix << "max_ghost_cm_error " << c.maxGhostCMError << "\n"
                        << prefix << "total_volume_error " << totalVolumeError << "\n";
                if(not ok)
                    metrics << "failed " << prefix << "comparison\n";
            }
        };
        const auto compare = [&](VoronoiGrid &grid, const std::string &name) { compareWith(grid, reference, name); };
        // A comparison that must hold bit for bit (up to the numbering of the ghosts, or not)
        const auto compareExactly = [&](VoronoiGrid &a, VoronoiGrid &b, const std::string &name, bool sameGhostNumbering)
        {
            const Differences local = sameGhostNumbering ? CountExactDifferences(a, b) : CountDifferencesUpToGhostNumbering(a, b);
            // (every rank has the same categories: a rank that stops at the sizes adds the rest as zeros)
            std::string what;
            std::size_t differences = 0;
            for(const char *category : {"cells", "mesh_points", "own_points", "faces", "face_vertices", "ghost_records", "sent_procs",
                                        "sent_points", "volumes", "centroids", "cell_faces", "face_neighbors", "face_areas", "face_centroids"})
            {
                const auto it = local.counts.find(category);
                const std::size_t global = GlobalSum(it == local.counts.end() ? 0 : it->second);
                differences += global;
                if(global > 0)
                    what += std::string(" ") + category + "=" + std::to_string(global);
            }
            pass = pass && (differences == 0);
            if(rank == 0)
            {
                std::cout << kind << " " << name << ": differences=" << differences << what << (differences == 0 ? "" : " FAIL") << std::endl;
                metrics << kind << "_" << name << "_differences " << differences << "\n";
                if(differences != 0)
                    metrics << "failed " << kind << "_" << name << "\n";
            }
        };

        for(std::size_t count : {1, 2, 4, 7})
        {
            VoronoiGrid grid(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
            grid.SetThreadSubdomains(count);
            BuildTwice(grid, points, kind);
            compare(grid, "subdomains=" + std::to_string(count));
            if(count == 1)
                compareExactly(grid, reference, "subdomains=1_vs_one_domain", false);
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

        // Some ranks ask for subdomains, the others do not: every rank builds as one domain
        if(worldSize > 1)
        {
            VoronoiGrid grid(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
            grid.SetThreadSubdomains((rank % 2 == 1) ? 4 : 0);
            BuildTwice(grid, points, kind);
            compareExactly(grid, reference, "mixed_requests_vs_one_domain", true);
        }

        // Load balance switches: with subdomains, MockMesh rebuilds from them with the ghosts
        // it knows; each step is compared with the same steps on a mesh built as one domain
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
                add(count, "subdomains=" + std::to_string(count));
            add(4, "subdomains=4_repeated");
            if(worldSize > 1)
                add((rank % 2 == 1) ? 4 : 0, "mixed_requests");
            for(std::size_t step = 0; step < SWITCH_STEPS; ++step)
            {
                for(std::unique_ptr<LoadBalanceSwitches> &s : switches)
                    s->Step(step, kind);
                VoronoiGrid &oneDomain = switches.front()->grid;
                const std::string at = std::string("switch_") + SWITCH_STEP_NAMES[step] + "_";
                if(step == 0)
                {
                    // the weights must move points between ranks
                    const std::size_t notRebalanced = GlobalSum(switches.front()->rebalanced ? 0 : 1);
                    pass = pass && (notRebalanced == 0);
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
                    else if(names[i] == "mixed_requests")
                        compareExactly(grid, oneDomain, at + "mixed_requests_vs_one_domain", true);
                    else
                        compareWith(grid, oneDomain, at + names[i]);
                }
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
