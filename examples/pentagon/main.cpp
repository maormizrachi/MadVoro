#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <mpi.h>
#include <Voronoi3D.hpp>
#include "../Vector3D.hpp"
#include "pentagon_helpers.hpp"

using Face = MadVoro::Face<Vector3D>;

#define SIDES 80

int main(int argc, char *argv[])
{
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    double R = 10;
    double theta = 2*M_PI / SIDES;
    double height = 3;

    std::vector<Vector3D> points_bottom;
    std::vector<Vector3D> points_top;

    for(size_t i = 0; i < SIDES; i++)
    {
        double x = cos(theta*i) * R;
        double y = sin(theta*i) * R;
        points_bottom.push_back(Vector3D(x, y, 0));
        points_top.push_back(Vector3D(x, y, height));
    }

    Face::point_vec_v base1_points, base2_points;
    for(const Vector3D &p : points_bottom)
    {
        base1_points.push_back(p);
    }
    for(const Vector3D &p : points_top)
    {
        base2_points.push_back(p);
    }
    Face base1(base1_points);
    Face base2(base2_points);
    
    std::vector<Face> box_faces = {base1, base2};
    for(size_t i = 0; i < SIDES; i++)
    {
        Face side({points_bottom[i], points_bottom[(i+1)%SIDES], points_top[(i+1)%SIDES], points_top[i]});
        box_faces.push_back(side);
    }

    MadVoro::Voronoi3D<Vector3D> diagram(box_faces);
    auto [ll, ur] = diagram.GetBoxCoordinates();

    if(rank == 0)
    {
        std::cout << "Generating points" << std::endl;
    }

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
    diagram.ToVTK("pentagon_example.vtu");
#else
    if(rank == 0)
    {
        std::cout << "VTK support not enabled. Skipping VTK output." << std::endl;
    }
#endif
    
    MPI_Finalize();

    return 0;
}
