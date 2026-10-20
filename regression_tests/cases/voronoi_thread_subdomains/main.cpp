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
#include <algorithm>
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

void BuildTwice(VoronoiGrid &grid, const std::vector<Vector3D> &points, const std::string &kind)
{
    MadVoro::regression_tests::BuildVoronoiMesh(grid, points);
    MadVoro::regression_tests::BuildVoronoiMesh(grid, Moved(points, kind));
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
    double totalVolume = 0;
    for(std::size_t c = 0; c < cells; ++c)
    {
        const double width = std::cbrt(reference.GetVolume(c));
        result.neighborMismatches += (Neighbors(subdomains, c, width) != Neighbors(reference, c, width)) ? 1 : 0;
        result.maxVolumeError = std::max(result.maxVolumeError,
                                         std::abs(subdomains.GetVolume(c) - reference.GetVolume(c)) / reference.GetVolume(c));
        result.maxCMError = std::max(result.maxCMError, Distance(subdomains.GetCellCM(c), reference.GetCellCM(c)) / width);
        totalVolume += subdomains.GetVolume(c);
    }
    result.totalVolumeError = MadVoro::regression_tests::ReduceGlobalVolume(totalVolume);
    result.totalVolumeError = std::abs(result.totalVolumeError - 1.0);
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
    for(const std::string kind : {"uniform", "lattice"})
    {
        const std::vector<Vector3D> points = MakePoints(kind, numPoints, rank);
        VoronoiGrid reference(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
        BuildTwice(reference, points, kind);

        VoronoiGrid one(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
        one.SetThreadSubdomains(1);
        BuildTwice(one, points, kind);
        const std::size_t exactDifferences = CountExactDifferences(one, reference);
        pass = pass && (exactDifferences == 0);
        if(rank == 0)
        {
            std::cout << kind << " subdomains=1: differences=" << exactDifferences << std::endl;
            metrics << kind << "_1_differences " << exactDifferences << "\n";
        }

        for(std::size_t count : {2, 4, 8})
        {
            VoronoiGrid grid(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
            grid.SetThreadSubdomains(count);
            BuildTwice(grid, points, kind);
            const Comparison c = Compare(grid, reference);
            const bool ok = c.neighborMismatches == 0 && c.unlistedFaces == 0 && c.maxVolumeError < MOMENTS_TOLERANCE &&
                            c.maxCMError < MOMENTS_TOLERANCE && c.totalVolumeError < MOMENTS_TOLERANCE;
            pass = pass && ok;
            if(rank == 0)
            {
                std::cout << kind << " subdomains=" << count << ": neighbor_mismatches=" << c.neighborMismatches
                          << " unlisted_faces=" << c.unlistedFaces << " max_volume_error=" << c.maxVolumeError
                          << " max_cm_error=" << c.maxCMError << " total_volume_error=" << c.totalVolumeError
                          << (ok ? "" : " FAIL") << std::endl;
                metrics << kind << "_" << count << "_neighbor_mismatches " << c.neighborMismatches << "\n"
                        << kind << "_" << count << "_unlisted_faces " << c.unlistedFaces << "\n"
                        << kind << "_" << count << "_max_volume_error " << c.maxVolumeError << "\n"
                        << kind << "_" << count << "_max_cm_error " << c.maxCMError << "\n"
                        << kind << "_" << count << "_total_volume_error " << c.totalVolumeError << "\n";
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
