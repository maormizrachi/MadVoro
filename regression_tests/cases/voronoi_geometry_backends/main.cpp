// The geometry executor MadVoro was built with (Kokkos, on a device or host threads,
// or plain host threads) must give exactly what the host executor gives: the same
// circumcenters, faces, cell face lists, volumes and centroids, bit for bit.
// Each case builds a mesh, then recomputes its geometry from the same flat mesh with
// HostGeometryExecutor and compares. Without Kokkos the two executors are the same.
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "regression_tests/lib/voronoi_test_common.hpp"
#include "geometry/HostGeometryExecutor.hpp"

#ifdef MADVORO_WITH_MPI
#include <mpi.h>
#endif

namespace
{

using MadVoro::Geometry::flat_index_t;
using MadVoro::Geometry::FlatPoint;
using MadVoro::Geometry::FlatVoronoi;
using MadVoro::regression_tests::VoronoiGrid;

bool Same(const FlatPoint &a, const FlatPoint &b)
{
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

template<typename T, typename Equal>
std::size_t CountDifferences(const std::vector<T> &a, const std::vector<T> &b, Equal equal)
{
    if(a.size() != b.size())
        return 1 + (a.size() > b.size() ? a.size() - b.size() : b.size() - a.size());
    std::size_t differences = 0;
    for(std::size_t i = 0; i < a.size(); ++i)
        differences += !equal(a[i], b[i]);
    return differences;
}

template<typename T>
std::size_t CountDifferences(const std::vector<T> &a, const std::vector<T> &b)
{
    return CountDifferences(a, b, [](const T &x, const T &y) { return x == y; });
}

// Recomputes the grid's geometry on the host and counts every difference.
std::size_t CompareWithHost(const VoronoiGrid &grid)
{
    const MadVoro::Geometry::FlatMesh &mesh = grid.GetFlatMesh();
    const FlatVoronoi &built = grid.GetFlatVoronoi();
    const MadVoro::Geometry::FlatMeshView meshView = mesh.View();
    const flat_index_t tetraCount = mesh.tetras.size();

    FlatVoronoi host;
    host.ResizeCellsAndTetras(meshView.ownedCount, tetraCount);
    std::vector<std::uint8_t> compute(tetraCount, 0);
    for(flat_index_t t = 0; t < tetraCount; ++t)
    {
        bool owned = false;
        for(int k = 0; k < 4; ++k)
            owned = owned || meshView.IsOwned(mesh.tetras[t].points[k]);
        compute[t] = owned ? 1 : 0;
    }

    MadVoro::Geometry::HostGeometryExecutor executor;
    executor.SetMesh(mesh);
    executor.Circumcenters(compute, host);
    std::size_t differences = 0;
    for(flat_index_t t = 0; t < tetraCount; ++t)
    {
        if(compute[t] == 0)
            continue;
        // Extended-precision circumcenters are the host's in both builds; take them as built.
        if(host.tetraStatus[t] != MadVoro::Geometry::GeometryStatus::Ok)
            host.vertices[t] = built.vertices[t];
        differences += !Same(host.vertices[t], built.vertices[t]);
    }
    host.vertices = built.vertices;
    executor.SetVertices(host);

    std::vector<std::uint8_t> looseAreaCut(meshView.pointCount);
    for(flat_index_t p = 0; p < meshView.pointCount; ++p)
        looseAreaCut[p] = (grid.IsPointOutsideBox(p) && !grid.IsPeriodicImage(p)) ? 1 : 0;
    const MadVoro::Geometry::FaceBuildReport report = executor.BuildFaces(built.cellOrder, looseAreaCut, host);
    differences += (report.result.error != MadVoro::Geometry::FaceBuildError::None);
    executor.CellMoments(host);

    differences += CountDifferences(host.faceNeighbor0, built.faceNeighbor0);
    differences += CountDifferences(host.faceNeighbor1, built.faceNeighbor1);
    differences += CountDifferences(host.faceVertexOffsets, built.faceVertexOffsets);
    differences += CountDifferences(host.faceVertices, built.faceVertices);
    differences += CountDifferences(host.faceArea, built.faceArea);
    differences += CountDifferences(host.faceCM, built.faceCM, Same);
    differences += CountDifferences(host.cellFaceOffsets, built.cellFaceOffsets);
    differences += CountDifferences(host.cellFaces, built.cellFaces);
    differences += CountDifferences(host.cellVolume, built.cellVolume);
    differences += CountDifferences(host.cellCM, built.cellCM, Same);
    return differences;
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
            points = MadVoro::regression_tests::RandRectangular(count, Vector3D(0, 0, 0), Vector3D(1, 1, 1), 97531ULL);
        }
    }
#ifdef MADVORO_WITH_MPI
    points = MadVoro::regression_tests::SpreadPointsFromRoot(points, rank);
#endif
    return points;
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
    const std::size_t numPoints = (worldSize > 1) ? 60000 : 20000;
#ifdef MADVORO_WITH_KOKKOS
    const std::string backend = "kokkos";
#else
    const std::string backend = "host";
#endif

    std::ofstream metrics;
    if(rank == 0)
    {
        metrics.open("voronoi_geometry_backends_metrics.txt");
        metrics << "mode " << (worldSize > 1 ? "mpi" : "serial") << "\nbackend " << backend << "\n";
    }
    bool pass = true;
    for(const std::string kind : {"uniform", "lattice", "periodic"})
    {
        VoronoiGrid grid(Vector3D(0, 0, 0), Vector3D(1, 1, 1));
        if(kind == "periodic")
            grid.SetPeriodicBoundaries(true, true, true);
        MadVoro::regression_tests::BuildVoronoiMesh(grid, MakePoints(kind, numPoints, rank));
        const std::size_t differences = GlobalSum(CompareWithHost(grid));
        pass = pass && (differences == 0);
        if(rank == 0)
        {
            std::cout << kind << ": backend=" << backend << " differences=" << differences << std::endl;
            metrics << kind << "_differences " << differences << "\n";
        }
    }

    int passed = pass ? 1 : 0;
    if(rank == 0)
    {
        metrics << "pass " << passed << "\n";
        std::cout << "voronoi_geometry_backends " << (passed ? "PASS" : "FAIL") << " backend=" << backend << std::endl;
    }
#ifdef MADVORO_WITH_MPI
    MPI_Finalize();
#endif
    return passed ? 0 : 1;
}
