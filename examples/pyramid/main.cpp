#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <mpi.h>
#include <Voronoi3D.hpp>
#include "../Vector3D.hpp"
#include "pyramid_helpers.hpp"

using Face = MadVoro::Face<Vector3D>;

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    double height = 15;

    Vector3D A(0, 0, 0), vec1(5, 0, 0), vec2(0, 5, 0);
    Vector3D B(A+vec1), C(A+vec1+vec2), D(A+vec2);

    Vector3D S = A + 0.5 * C + Vector3D(0, 0, height);

    Face base({A, B, C, D});
    Face side1({A, D, S});
    Face side2({C, D, S});
    Face side3({B, C, S});
    Face side4({A, B, S});
    std::vector<Face> box_faces = {base, side1, side2, side3, side4};

    MadVoro::Voronoi3D<Vector3D> diagram(box_faces);
    auto [ll, ur] = diagram.GetBoxCoordinates();

    size_t N = 1000;
    std::vector<Vector3D> points = GeneratePoints(box_faces, N, ll, ur);

    if(rank == 0)
    {
        std::cout << "Starting voronoi build" << std::endl;
    }
    points = diagram.BuildParallel(points);

    MPI_Barrier(MPI_COMM_WORLD);

#ifdef MADVORO_WITH_VTK
    if(rank == 0)
    {
        std::cout << "Starting VTK output" << std::endl;
    }
    diagram.ToVTK("pyramid_example.vtu");
#else
    if(rank == 0)
    {
        std::cout << "VTK support not enabled. Skipping VTK output." << std::endl;
    }
#endif
    
    MPI_Finalize();

    return 0;
}
