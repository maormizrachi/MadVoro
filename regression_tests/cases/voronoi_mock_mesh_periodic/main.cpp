#include <cmath>
#include <cstddef>
#include <fstream>
#include <iostream>
#include <vector>

#include "Voronoi3D.hpp"
#include "examples/Vector3D.hpp"

#include <mpi.h>

namespace
{

double SumOwnedVolume(const MadVoro::Voronoi3D<Vector3D> &voronoi)
{
    double localVolume = 0.0;
    for(std::size_t i = 0; i < voronoi.GetPointNo(); ++i)
    {
        localVolume += voronoi.GetVolume(i);
    }
    double totalVolume = localVolume;
    MPI_Allreduce(&localVolume, &totalVolume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return totalVolume;
}

std::vector<Vector3D> BuildPeriodicLatticePoints(int rank, int worldSize)
{
    std::vector<Vector3D> localPoints;
    for(int ix = 0; ix < 2; ++ix)
    {
        for(int iy = 0; iy < 2; ++iy)
        {
            for(int iz = 0; iz < 2; ++iz)
            {
                std::size_t globalIdx = static_cast<std::size_t>(ix * 4 + iy * 2 + iz);
                if(globalIdx % static_cast<std::size_t>(worldSize) == static_cast<std::size_t>(rank))
                {
                    localPoints.emplace_back(0.25 + 0.5 * ix, 0.25 + 0.5 * iy, 0.25 + 0.5 * iz);
                }
            }
        }
    }
    return localPoints;
}

} // namespace

int main()
{
    int rank = 0;
    int worldSize = 1;
    MPI_Init(nullptr, nullptr);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);

    const Vector3D ll(0.0, 0.0, 0.0);
    const Vector3D ur(1.0, 1.0, 1.0);
    const double boxVolume = (ur.x - ll.x) * (ur.y - ll.y) * (ur.z - ll.z);

    MadVoro::Voronoi3D<Vector3D> voronoi(ll, ur);
    voronoi.SetPeriodicBoundaries(true, true, true);
    voronoi.BuildParallel(BuildPeriodicLatticePoints(rank, worldSize));

    const double volumeBefore = SumOwnedVolume(voronoi);

    std::vector<double> weights(voronoi.getAllPoints().size(), 1.0);
    for(std::size_t i = 0; i < weights.size(); ++i)
    {
        weights[i] = static_cast<double>((i * 17 + static_cast<std::size_t>(rank) * 11) % 23 + 1);
    }
    voronoi.Rebalance(weights);

    const double volumeAfter = SumOwnedVolume(voronoi);
    const double relErrorBefore = std::abs(volumeBefore - boxVolume) / boxVolume;
    const double relErrorAfter = std::abs(volumeAfter - boxVolume) / boxVolume;
    const int mockMeshUsed = voronoi.DidRebalance() ? 1 : 0;
    const int passed = (relErrorBefore < 1e-10 && relErrorAfter < 1e-10 && mockMeshUsed == 1) ? 1 : 0;

    if(rank == 0)
    {
        std::cout << "voronoi_mock_mesh_periodic PASS=" << (passed ? 1 : 0)
                  << " volume_before=" << volumeBefore
                  << " volume_after=" << volumeAfter
                  << " rel_error_before=" << relErrorBefore
                  << " rel_error_after=" << relErrorAfter
                  << " did_rebalance=" << mockMeshUsed
                  << std::endl;

        std::ofstream out("voronoi_mock_mesh_periodic_metrics.txt");
        out.setf(std::ios::scientific);
        out.precision(16);
        out << "mode mpi\n";
        out << "volume_before " << volumeBefore << "\n";
        out << "volume_after " << volumeAfter << "\n";
        out << "rel_error_before " << relErrorBefore << "\n";
        out << "rel_error_after " << relErrorAfter << "\n";
        out << "did_rebalance " << mockMeshUsed << "\n";
        out << "pass " << passed << "\n";
    }

    MPI_Finalize();
    return passed ? 0 : 1;
}
