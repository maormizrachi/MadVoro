// Building from thread subdomains must give the tessellation of a build as one domain.
// With one subdomain the result must be the same bit for bit (points, faces, face
// vertices, volumes, centroids, radii carried to the next build). With several, each
// cell must have the same neighbors and its volume and centroid must agree closely.
// Every face between two cells must be listed by both. Each case builds twice: on the
// points, then on slightly moved points, so the second build starts from the radii the
// first one left.
// Near a wall a cell, its mirror, another point and that point's mirror are always on
// one sphere, so a build may leave a sliver face (of zero area in exact arithmetic)
// between a cell and the mirror of another point. A subdomain only mirrors its own
// points, so its slivers differ from those of a build as one domain: faces against
// mirrors are compared above SLIVER_AREA, and volumes and centroids to MOMENTS_TOLERANCE,
// orders of magnitude below the error of a missing neighbor.
// The same comparison holds in a periodic box (all axes, or x only) and for partial
// builds (a scattered third of the points, also in a periodic box), where a subdomain
// adds its ghosts in another order than one domain does, so one subdomain is compared
// within the tolerances too; there, building again with 4 subdomains must give the same
// result bit for bit.
#include <algorithm>
#include <array>
#include <cstdint>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <set>
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

// Relative to the cell's width squared
constexpr double CELL_FACE_AREA = 1e-8;
constexpr double SLIVER_AREA = 1e-6;
constexpr double MOMENTS_TOLERANCE = 1e-7;

bool Same(const Vector3D &a, const Vector3D &b)
{
    return a.x == b.x && a.y == b.y && a.z == b.z;
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
        else
        {
            points = MadVoro::regression_tests::RandRectangular(count, Vector3D(0, 0, 0), Vector3D(1, 1, 1), 24680ULL);
        }
    }
#ifdef MADVORO_WITH_MPI
    points = MadVoro::regression_tests::SpreadPointsFromRoot(points, rank);
#endif
    return points;
}

// Small deterministic motion (the lattice stays fixed)
std::vector<Vector3D> Moved(const std::vector<Vector3D> &points, const std::string &kind)
{
    if(kind == "lattice")
        return points;
    std::vector<Vector3D> moved = points;
    for(std::size_t i = 0; i < moved.size(); ++i)
    {
        const double phase = static_cast<double>(i);
        moved[i].x = std::min(std::max(moved[i].x + 1e-3 * std::sin(phase), 1e-9), 1.0 - 1e-9);
        moved[i].y = std::min(std::max(moved[i].y + 1e-3 * std::cos(1.3 * phase), 1e-9), 1.0 - 1e-9);
        moved[i].z = std::min(std::max(moved[i].z + 1e-3 * std::sin(0.7 * phase + 1.0), 1e-9), 1.0 - 1e-9);
    }
    return moved;
}

// The box's boundaries and which points are built
struct Setup
{
    std::string name;
    std::array<bool, 3> periodic;
    bool partial;
};

// A scattered third of the points, by position
bool InScatteredThird(const Vector3D &point)
{
    const std::uint64_t qx = static_cast<std::uint64_t>(point.x * 1048576.0);
    const std::uint64_t qy = static_cast<std::uint64_t>(point.y * 1048576.0);
    const std::uint64_t qz = static_cast<std::uint64_t>(point.z * 1048576.0);
    return ((qx * 73856093ULL) ^ (qy * 19349663ULL) ^ (qz * 83492791ULL)) % 3 == 0;
}

void Build(VoronoiGrid &grid, const std::vector<Vector3D> &points, const std::vector<std::size_t> &active, const Setup &setup)
{
    if(not setup.partial)
    {
        MadVoro::regression_tests::BuildVoronoiMesh(grid, points);
        return;
    }
#ifdef MADVORO_WITH_MPI
    grid.BuildPartiallyParallel(points, active);
#else
    grid.BuildPartially(points, active);
#endif
}

void BuildTwice(VoronoiGrid &grid, const std::vector<Vector3D> &points, const std::string &kind, const Setup &setup)
{
    grid.SetPeriodicBoundaries(setup.periodic[0], setup.periodic[1], setup.periodic[2]);
    std::vector<std::size_t> active;
    for(std::size_t i = 0; i < points.size(); ++i)
    {
        if(InScatteredThird(points[i]))
            active.push_back(i);
    }
    Build(grid, points, active, setup);
    Build(grid, Moved(points, kind), active, setup);
}

// Differences between two builds that must be the same bit for bit
std::size_t CountExactDifferences(VoronoiGrid &a, VoronoiGrid &b)
{
    std::size_t differences = 0;
    const auto count = [&differences](bool different) { differences += different ? 1 : 0; };
    count(a.GetPointNo() != b.GetPointNo());
    count(a.getMeshPoints().size() != b.getMeshPoints().size());
    count(a.GetTotalFacesNumber() != b.GetTotalFacesNumber());
    count(a.GetFacePoints().size() != b.GetFacePoints().size());
    if(differences > 0)
        return differences;
    for(std::size_t p = 0; p < a.getMeshPoints().size(); ++p)
        count(!Same(a.getMeshPoints()[p], b.getMeshPoints()[p]));
    for(std::size_t c = 0; c < a.GetPointNo(); ++c)
    {
        count(a.GetVolume(c) != b.GetVolume(c));
        count(!Same(a.GetCellCM(c), b.GetCellCM(c)));
        const auto &facesA = a.GetCellFaces(c);
        const auto &facesB = b.GetCellFaces(c);
        count(!std::equal(facesA.begin(), facesA.end(), facesB.begin(), facesB.end()));
    }
    for(std::size_t f = 0; f < a.GetTotalFacesNumber(); ++f)
    {
        count(a.GetFaceNeighbors(f) != b.GetFaceNeighbors(f));
        count(a.GetArea(f) != b.GetArea(f));
        count(!Same(a.FaceCM(f), b.FaceCM(f)));
        const auto &verticesA = a.GetPointsInFace(f);
        const auto &verticesB = b.GetPointsInFace(f);
        count(verticesA.size() != verticesB.size());
        for(std::size_t k = 0; k < std::min(verticesA.size(), verticesB.size()); ++k)
            count(!Same(a.GetFacePoints()[verticesA[k]], b.GetFacePoints()[verticesB[k]]));
    }
    for(std::size_t p = 0; p < a.getMeshPoints().size(); ++p)
        count(!Same(a.GetAllCM()[p], b.GetAllCM()[p]));
    return differences;
}

// A cell's neighbors over faces that are not negligibly small: own cells by index,
// ghosts by position
std::set<std::tuple<std::size_t, double, double, double>> Neighbors(const VoronoiGrid &grid, std::size_t cell, double width)
{
    std::set<std::tuple<std::size_t, double, double, double>> neighbors;
    for(std::size_t face : grid.GetCellFaces(cell))
    {
        const auto &pair = grid.GetFaceNeighbors(face);
        const std::size_t other = (pair.first == cell) ? pair.second : pair.first;
        const double minArea = ((other < grid.GetPointNo()) ? CELL_FACE_AREA : SLIVER_AREA) * width * width;
        if(grid.GetArea(face) <= minArea)
            continue;
        if(other < grid.GetPointNo())
        {
            neighbors.emplace(other, 0.0, 0.0, 0.0);
        }
        else
        {
            const Vector3D &p = grid.GetMeshPoint(other);
            neighbors.emplace(grid.GetPointNo(), p.x, p.y, p.z);
        }
    }
    return neighbors;
}

struct Comparison
{
    std::size_t neighborMismatches = 0;
    std::size_t unlistedFaces = 0;
    double maxVolumeError = 0;
    double maxCMError = 0;
    double totalVolumeError = 0;
};

// Compares a subdomain build with a build as one domain
Comparison Compare(const VoronoiGrid &subdomains, const VoronoiGrid &reference)
{
    Comparison result;
    const std::size_t cells = reference.GetPointNo();
    if(subdomains.GetPointNo() != cells)
    {
        result.neighborMismatches = cells + 1;
        return result;
    }
    double totalVolume = 0, referenceVolume = 0;
    for(std::size_t c = 0; c < cells; ++c)
    {
        const double width = std::cbrt(reference.GetVolume(c));
        result.neighborMismatches += (Neighbors(subdomains, c, width) != Neighbors(reference, c, width)) ? 1 : 0;
        result.maxVolumeError = std::max(result.maxVolumeError,
                                         std::abs(subdomains.GetVolume(c) - reference.GetVolume(c)) / reference.GetVolume(c));
        result.maxCMError = std::max(result.maxCMError, Distance(subdomains.GetCellCM(c), reference.GetCellCM(c)) / width);
        totalVolume += subdomains.GetVolume(c);
        referenceVolume += reference.GetVolume(c);
    }
    // (the whole box in a full build)
    referenceVolume = MadVoro::regression_tests::ReduceGlobalVolume(referenceVolume);
    result.totalVolumeError = std::abs(MadVoro::regression_tests::ReduceGlobalVolume(totalVolume) - referenceVolume) / referenceVolume;
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
    return result;
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
    const std::size_t numPoints = 20000;

    std::ofstream metrics;
    if(rank == 0)
    {
        metrics.open("voronoi_thread_subdomains_metrics.txt");
        metrics << "mode " << (worldSize > 1 ? "mpi" : "serial") << "\n";
    }
    bool pass = (worldSize == 1);
    if(rank == 0 && worldSize != 1)
        std::cout << "voronoi_thread_subdomains runs on one rank" << std::endl;
    const std::vector<Setup> setups = {{"box", {false, false, false}, false},
                                       {"periodic", {true, true, true}, false},
                                       {"periodic_x", {true, false, false}, false},
                                       {"partial", {false, false, false}, true},
                                       {"partial_periodic", {true, true, true}, true}};
    for(const Setup &setup : setups)
    for(const std::string kind : {"uniform", "lattice"})
    {
        // (the box keeps the names of its metrics)
        const std::string name = (setup.name == "box") ? kind : setup.name + "_" + kind;
        const std::vector<Vector3D> points = MakePoints(kind, numPoints, rank);
        VoronoiGrid reference(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
        BuildTwice(reference, points, kind, setup);

        if(setup.name == "box")
        {
            VoronoiGrid one(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
            one.SetThreadSubdomains(1);
            BuildTwice(one, points, kind, setup);
            const std::size_t exactDifferences = CountExactDifferences(one, reference);
            pass = pass && (exactDifferences == 0);
            if(rank == 0)
            {
                std::cout << name << " subdomains=1: differences=" << exactDifferences << std::endl;
                metrics << name << "_1_differences " << exactDifferences << "\n";
            }
        }

        const std::vector<std::size_t> counts = (setup.name == "box") ? std::vector<std::size_t>{2, 4, 8} : std::vector<std::size_t>{1, 2, 4, 8};
        for(std::size_t count : counts)
        {
            VoronoiGrid grid(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
            grid.SetThreadSubdomains(count);
            BuildTwice(grid, points, kind, setup);
            const Comparison c = Compare(grid, reference);
            const bool ok = c.neighborMismatches == 0 && c.unlistedFaces == 0 && c.maxVolumeError < MOMENTS_TOLERANCE &&
                            c.maxCMError < MOMENTS_TOLERANCE && c.totalVolumeError < MOMENTS_TOLERANCE;
            pass = pass && ok;
            if(rank == 0)
            {
                std::cout << name << " subdomains=" << count << ": neighbor_mismatches=" << c.neighborMismatches
                          << " unlisted_faces=" << c.unlistedFaces << " max_volume_error=" << c.maxVolumeError
                          << " max_cm_error=" << c.maxCMError << " total_volume_error=" << c.totalVolumeError
                          << (ok ? "" : " FAIL") << std::endl;
                metrics << name << "_" << count << "_neighbor_mismatches " << c.neighborMismatches << "\n"
                        << name << "_" << count << "_unlisted_faces " << c.unlistedFaces << "\n"
                        << name << "_" << count << "_max_volume_error " << c.maxVolumeError << "\n"
                        << name << "_" << count << "_max_cm_error " << c.maxCMError << "\n"
                        << name << "_" << count << "_total_volume_error " << c.totalVolumeError << "\n";
            }
            if(count == 4 and setup.name != "box")
            {
                VoronoiGrid again(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
                again.SetThreadSubdomains(count);
                BuildTwice(again, points, kind, setup);
                const std::size_t repeatedDifferences = CountExactDifferences(again, grid);
                pass = pass && (repeatedDifferences == 0);
                if(rank == 0)
                {
                    std::cout << name << " subdomains=4_repeated: differences=" << repeatedDifferences
                              << (repeatedDifferences == 0 ? "" : " FAIL") << std::endl;
                    metrics << name << "_4_repeated_differences " << repeatedDifferences << "\n";
                }
            }
        }
    }

    const int passed = pass ? 1 : 0;
    if(rank == 0)
    {
        metrics << "pass " << passed << "\n";
        std::cout << "voronoi_thread_subdomains " << (passed ? "PASS" : "FAIL") << std::endl;
    }
#ifdef MADVORO_WITH_MPI
    MPI_Finalize();
#endif
    return passed ? 0 : 1;
}
