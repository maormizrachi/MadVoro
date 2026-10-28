// Checks that FlatMesh and FlatVoronoi describe exactly the Delaunay3D / Voronoi3D they were built from.
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "regression_tests/lib/voronoi_test_common.hpp"
#include "geometry/FlatMesh.hpp"
#include "geometry/FlatVoronoi.hpp"

#ifdef MADVORO_WITH_MPI
#include <mpi.h>
#endif

namespace
{

using MadVoro::Geometry::flat_index_t;
using MadVoro::regression_tests::VoronoiGrid;

bool Same(const MadVoro::Geometry::FlatPoint &a, const Vector3D &b)
{
    return a.x == b.x and a.y == b.y and a.z == b.z;
}

// Returns the number of violations found in the flat Delaunay mesh.
std::size_t CheckFlatMesh(const MadVoro::Delaunay3D<Vector3D> &del)
{
    std::size_t errors = 0;
    const MadVoro::Geometry::FlatMesh mesh = MadVoro::Geometry::MakeFlatMesh(del);
    const MadVoro::Geometry::FlatMeshView view = mesh.View();

    errors += (view.pointCount != del.points_.size()) + (view.tetraCount != del.tetras_.size()) + (view.ownedCount != del.Norg_);
    for(flat_index_t i = 0; i < view.pointCount; ++i)
    {
        errors += not Same(view.points[i], del.points_[i]);
    }
    for(const std::size_t &t : del.empty_tetras_)
    {
        errors += view.IsAlive(t);
    }

    std::size_t expectedListEntries = 0;
    for(flat_index_t t = 0; t < view.tetraCount; ++t)
    {
        if(not view.IsAlive(t))
        {
            continue;
        }
        for(int k = 0; k < 4; ++k)
        {
            const flat_index_t p = view.tetras[t].points[k];
            errors += (p >= view.pointCount);
            if(view.IsOwned(p))
            {
                ++expectedListEntries;
                // t must appear exactly once in p's list
                std::size_t found = 0;
                for(flat_index_t j = 0; j < view.PointTetraCount(p); ++j)
                {
                    found += (view.PointTetra(p, j) == t);
                }
                errors += (found != 1);
            }
            // a neighbour inside the mesh must point back at t
            const flat_index_t n = view.tetras[t].neighbors[k];
            if(n < view.tetraCount and view.IsAlive(n))
            {
                const flat_index_t *back = view.tetras[n].neighbors;
                errors += not (back[0] == t or back[1] == t or back[2] == t or back[3] == t);
            }
        }
    }
    errors += (view.pointTetraOffsets[0] != 0) + (view.pointTetraOffsets[view.ownedCount] != expectedListEntries);
    for(flat_index_t p = 0; p < view.ownedCount; ++p)
    {
        errors += (view.pointTetraOffsets[p + 1] < view.pointTetraOffsets[p]);
        for(flat_index_t j = 0; j < view.PointTetraCount(p); ++j)
        {
            const flat_index_t t = view.PointTetra(p, j);
            const flat_index_t *v = view.tetras[t].points;
            errors += not (view.IsAlive(t) and (v[0] == p or v[1] == p or v[2] == p or v[3] == p));
        }
    }

    // The overload taking lists must keep them as given, order included.
    std::vector<std::vector<std::size_t>> reversed(view.ownedCount);
    for(flat_index_t p = 0; p < view.ownedCount; ++p)
    {
        for(flat_index_t j = view.PointTetraCount(p); j-- > 0;)
        {
            reversed[p].push_back(view.PointTetra(p, j));
        }
    }
    const MadVoro::Geometry::FlatMesh kept = MadVoro::Geometry::MakeFlatMesh(del, reversed);
    for(flat_index_t p = 0; p < view.ownedCount; ++p)
    {
        const MadVoro::Geometry::FlatMeshView keptView = kept.View();
        errors += (keptView.PointTetraCount(p) != reversed[p].size());
        for(flat_index_t j = 0; j < keptView.PointTetraCount(p); ++j)
        {
            errors += (keptView.PointTetra(p, j) != reversed[p][j]);
        }
    }
    return errors;
}

// Returns the number of violations found in the flat Voronoi tessellation.
std::size_t CheckFlatVoronoi(const VoronoiGrid &grid, double &volumeSum)
{
    std::size_t errors = 0;
    const MadVoro::Geometry::FlatVoronoi flat = MadVoro::Geometry::MakeFlatVoronoi(grid);
    const MadVoro::Geometry::ConstFlatVoronoiView view = flat.View();

    errors += (view.cellCount != grid.GetPointNo()) + (view.faceCount != grid.GetTotalFacesNumber()) + (view.tetraCount != grid.GetFacePoints().size());
    for(flat_index_t t = 0; t < view.tetraCount; ++t)
    {
        errors += not Same(view.vertices[t], grid.GetFacePoints()[t]);
    }

    volumeSum = 0;
    errors += (view.cellFaceOffsets[0] != 0) + (view.faceVertexOffsets[0] != 0);
    for(flat_index_t c = 0; c < view.cellCount; ++c)
    {
        const auto &faces = grid.GetCellFaces(c);
        errors += (view.CellFaceCount(c) != faces.size());
        for(flat_index_t k = 0; k < view.CellFaceCount(c) and k < faces.size(); ++k)
        {
            const flat_index_t f = view.CellFace(c, k);
            errors += (f != faces[k]);
            // the face must separate c from the neighbour OtherSide reports
            errors += not (view.faceNeighbor0[f] == c or view.faceNeighbor1[f] == c);
            const auto &n = grid.GetFaceNeighbors(f);
            errors += (view.OtherSide(f, c) != (n.first == c ? n.second : n.first));
        }
        errors += (view.cellVolume[c] != grid.GetVolume(c)) + not Same(view.cellCM[c], grid.GetCellCM(c));
        volumeSum += view.cellVolume[c];
    }
    for(flat_index_t f = 0; f < view.faceCount; ++f)
    {
        const auto &vertices = grid.GetPointsInFace(f);
        errors += (view.FaceVertexCount(f) != vertices.size());
        for(flat_index_t k = 0; k < view.FaceVertexCount(f) and k < vertices.size(); ++k)
        {
            errors += (view.FaceVertex(f, k) != vertices[k]) + (view.FaceVertex(f, k) >= view.tetraCount);
        }
        const auto &n = grid.GetFaceNeighbors(f);
        errors += (view.faceNeighbor0[f] != n.first) + (view.faceNeighbor1[f] != n.second);
        errors += (view.faceArea[f] != grid.GetArea(f)) + not Same(view.faceCM[f], grid.FaceCM(f));
    }
    for(flat_index_t c = 0; c < view.cellCount; ++c)
    {
        errors += (view.cellStatus[c] != MadVoro::Geometry::GeometryStatus::Ok);
    }
    return errors;
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
    const std::uint64_t seed = 246813579ULL;
    const std::size_t numPoints = (worldSize > 1) ? static_cast<std::size_t>(40000) : static_cast<std::size_t>(10000);

    // Flat Delaunay mesh, on every rank from its own points (no ghosts needed).
    std::size_t meshErrors = 0;
    {
        std::vector<Vector3D> points = MadVoro::regression_tests::RandRectangular(4000, ll, ur, seed + static_cast<std::uint64_t>(rank));
        std::vector<std::size_t> order = HilbertOrder3D(points);
        MadVoro::Delaunay3D<Vector3D> del;
        del.Build(points, ur, ll, order);
        meshErrors = GlobalSum(CheckFlatMesh(del));
    }

    // Flat Voronoi, from a full (parallel, under MPI) build.
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
    VoronoiGrid grid(ll, ur);
    MadVoro::regression_tests::BuildVoronoiMesh(grid, points);
    double localVolume = 0;
    const std::size_t voronoiErrors = GlobalSum(CheckFlatVoronoi(grid, localVolume));
    const double totalVolume = MadVoro::regression_tests::ReduceGlobalVolume(localVolume);
    const double relError = std::abs(totalVolume - 1.0);

    const int passed = (meshErrors == 0 and voronoiErrors == 0 and relError < 1e-10) ? 1 : 0;
    if(rank == 0)
    {
        std::cout << "flat_mesh_errors=" << meshErrors << " flat_voronoi_errors=" << voronoiErrors << " volume_rel_error=" << relError << std::endl;
        std::ofstream metrics("voronoi_flat_layout_metrics.txt");
        metrics << "mode " << (worldSize > 1 ? "mpi" : "serial") << "\n" << "flat_mesh_errors " << meshErrors << "\n" << "flat_voronoi_errors " << voronoiErrors << "\n" << "volume_rel_error " << relError << "\n" << "pass " << passed << "\n";
        std::cout << "voronoi_flat_layout " << (passed ? "PASS" : "FAIL") << " seed=" << seed << std::endl;
    }

#ifdef MADVORO_WITH_MPI
    MPI_Finalize();
#endif
    return passed ? 0 : 1;
}
