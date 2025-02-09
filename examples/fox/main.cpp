#include <vector>
#include <mpi.h>
#include <algorithm>
#include <tuple> // for std::tie
#include <madvoro/Voronoi3D.hpp>
#include "mpi_utils.hpp"
#include "read_points.hpp"

using namespace MadVoro;

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // some of the points are inside the fox, and some of them are outside. This is crucial, of course, to create the fox's faces.
    // we should remember what points are inside.
    std::vector<Vector3D> allPoints;
    std::vector<double> isInside;
    
    if(rank == 0)
    {
        // read points
        std::vector<Vector3D> insidePoints = readPoints("input/points_inside_fox");
        std::vector<Vector3D> outsidePoints = readPoints("input/points_outside_fox");
        allPoints.insert(allPoints.end(), insidePoints.cbegin(), insidePoints.cend());
        for(size_t i = 0; i < insidePoints.size(); i++) isInside.push_back(1);
        allPoints.insert(allPoints.end(), outsidePoints.cbegin(), outsidePoints.cend());
        for(size_t i = 0; i < outsidePoints.size(); i++) isInside.push_back(0);
    }

    // not really necessary, but it is more efficient to distribute the points between the different ranks
    auto [myPoints, myIsInside] = SpreadPointsToProcessors(allPoints, isInside);

    Vector3D ll(-13, -1, -90), ur(15, 80, 70);
    Voronoi3D diag(ll, ur);
    
    if(rank == 0)
    {
        std::cout << "Constructing the Voronoi diagram" << std::endl;
    }
    diag.SetVerbosity(true); // print build iterations

    diag.BuildParallel(myPoints);

    // after the build, the list of points changes. The `myIsInside` list is not longer accurate, and should be updated.
    // we use a special function of Voronoi3D to help us recognize which points have moved, and where to.
    std::tie(myPoints, myIsInside) = GetPointsAfterBuildExchange(diag, myPoints, myIsInside);

    std::string vtkFilename = "fox.vtu";
    if(rank == 0)
    {
        std::cout << "Starting VTK Print to file called " << vtkFilename << std::endl;
    }

    // print to VTK, with a field of "isInside" to each point. When you want to see the VTK result, make sure you filter by 'isInside == 1'
    diag.ToVTK(vtkFilename, {"isInside"}, {myIsInside});

    MPI_Finalize();
}