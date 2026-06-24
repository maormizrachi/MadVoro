#include <vector>
#include <mpi.h>
#include <algorithm>
#include <tuple>
#include <madvoro/Voronoi3D.hpp>
#include "Vector3D.hpp"
#include "mpi_utils.hpp"
#include "read_points.hpp"

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::vector<Vector3D> allPoints;
    std::vector<double> isInside;
    
    if(rank == 0)
    {
        std::vector<Vector3D> insidePoints = readPoints("input/points_inside_fox");
        std::vector<Vector3D> outsidePoints = readPoints("input/points_outside_fox");
        allPoints.insert(allPoints.end(), insidePoints.cbegin(), insidePoints.cend());
        for(size_t i = 0; i < insidePoints.size(); i++) isInside.push_back(1);
        allPoints.insert(allPoints.end(), outsidePoints.cbegin(), outsidePoints.cend());
        for(size_t i = 0; i < outsidePoints.size(); i++) isInside.push_back(0);
    }

    auto [myPoints, myIsInside] = SpreadPointsToProcessors(allPoints, isInside);

    Vector3D ll(-13, -1, -90), ur(15, 80, 70);
    MadVoro::Voronoi3D<Vector3D> diag(ll, ur);
    
    if(rank == 0)
    {
        std::cout << "Constructing the Voronoi diagram" << std::endl;
    }
    diag.SetVerbosity(true);

    diag.BuildParallel(myPoints);

    std::tie(myPoints, myIsInside) = GetPointsAfterBuildExchange(diag, myPoints, myIsInside);

#ifdef MADVORO_WITH_VTK
    std::string vtkFilename = "fox.vtu";
    if(rank == 0)
    {
        std::cout << "Starting VTK Print to file called " << vtkFilename << std::endl;
    }
    diag.ToVTK(vtkFilename, {"isInside"}, {myIsInside});
#else
    if(rank == 0)
    {
        std::cout << "VTK support not enabled. Skipping VTK output." << std::endl;
    }
#endif

    MPI_Finalize();
}
