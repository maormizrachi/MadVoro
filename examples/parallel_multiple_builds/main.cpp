#include <iostream>
#include <chrono>
#include <mpi.h>
#include <Voronoi3D.hpp>
#include "../Vector3D.hpp"
#include "read_points.hpp"

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    std::string inputFileName = "input/" + std::to_string(rank);
    std::vector<Vector3D> myPoints = readPoints(inputFileName);
    
    const size_t iterations = 15;

    Vector3D ll(0, 0, 0), ur(1, 1, 1);
    MadVoro::Voronoi3D<Vector3D> distriburedDiagram(ll, ur);

    std::chrono::high_resolution_clock::time_point start, end;
    
    double totalTime = 0, firstTime = 0;
    for(size_t i = 0; i < iterations; ++i)
    {
        if(rank == 0)
        {
            std::cout << "Build " << i << ", ";
        }
        start = std::chrono::high_resolution_clock::now();
        myPoints = distriburedDiagram.BuildParallel(myPoints);
        end = std::chrono::high_resolution_clock::now();

        double time = std::chrono::duration<double>(end - start).count();
        if(rank == 0)
        {
            std::cout << time << " seconds" << std::endl;
        }
        if(i == 0) firstTime = time;
        if(i > 0) totalTime += time;
    }
    
    if(rank == 0)
    {
        std::cout << "Done with all builds!" << std::endl;
        std::cout << "First build took " << firstTime << " Later " << (iterations - 1) << " advanced builds took in average " << (totalTime) / (iterations - 1) << " seconds." << std::endl;
    }

    MPI_Finalize();
    return 0;
}
