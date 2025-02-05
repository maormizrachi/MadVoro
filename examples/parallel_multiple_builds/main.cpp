#include <iostream>
#include <chrono>
#include <mpi.h>
#include <madvoro/Voronoi3D.hpp>
#include "read_points.hpp"

using namespace MadVoro;

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // a list of all points
    std::string inputFileName = "input/" + std::to_string(rank);
    std::vector<Vector3D> myPoints = readPoints(inputFileName);

    Vector3D ll(0, 0, 0), ur(1, 1, 1);
    Voronoi3D distriburedDiagram(ll, ur);

    std::chrono::high_resolution_clock::time_point start, end;
    
    start = std::chrono::high_resolution_clock::now();
    // we re-assign to `myPoints` to get our new points list (after redistribution)
    myPoints = distriburedDiagram.BuildParallel(myPoints);
    end = std::chrono::high_resolution_clock::now();
    
    double time = std::chrono::duration<double>(end - start).count();
    
    if(rank == 0)
    {
        std::cout << "Done with building!" << std::endl;
        std::cout << "Single build time with " << size << " processes, is " << time << std::endl;
        std::cout << "Now writing to VTK" << std::endl;
    }

    distriburedDiagram.ToVTK("output.vtk");

    MPI_Finalize();
    return 0;
}