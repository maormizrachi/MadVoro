#include <iostream>
#include <vector>
#include <chrono>
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <mpi.h>
#include <madvoro/Voronoi3D.hpp>

#define DEFAULT_N 10000

using namespace MadVoro;

std::vector<Vector3D> GenerateRandomPoints(size_t N, const Vector3D &ll, const Vector3D &ur)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<Vector3D> result;

    boost::random::uniform_real_distribution<double> distribution_x(ll.x, ur.x);
    boost::random::uniform_real_distribution<double> distribution_y(ll.y, ur.y);
    boost::random::uniform_real_distribution<double> distribution_z(ll.z, ur.z);

    boost::random::mt19937 rng;
    rng.seed(rank);

    for(size_t i = 0; i < N; ++i)
    {
        result.push_back(Vector3D(distribution_x(rng), distribution_y(rng), distribution_z(rng)));
    }

    return result;
}

int main(int argc, char *argv[])
{
    if(argc != 1 and argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " [N]" << std::endl;
        return EXIT_FAILURE;
    }

    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    size_t N = (argc == 2)? std::stoul(argv[1]) : DEFAULT_N;

    Vector3D ll(-1, -1, -1), ur(1, 1, 1);
    std::vector<Vector3D> points = GenerateRandomPoints(N, ll, ur);
    
    MPI_Allreduce(MPI_IN_PLACE, &N, 1, MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);

    if(rank == 0)
    {
        std::cout << "In total, generated " << N << " random points (average per processor: " << N / size << ")" << std::endl;
        std::cout << "Starting build" << std::endl;
    }

    Voronoi3D voronoi(ll, ur);

    std::chrono::time_point<std::chrono::system_clock> start, end;

    start = std::chrono::system_clock::now();
    points = voronoi.BuildParallel(points);
    end = std::chrono::system_clock::now();
    
    double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    
    if(rank == 0)
    {
        std::cout << "End of build, time taken for first build is " << elapsed << " seconds" << std::endl;
        std::cout << "Rate (points/second, for a single processor): " << (points.size() / elapsed) << std::endl;
    }

    MPI_Finalize();
    return EXIT_SUCCESS;
}