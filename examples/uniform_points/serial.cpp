#include <iostream>
#include <vector>
#include <chrono>
#include <boost/random.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <madvoro/Voronoi3D.hpp>

#define DEFAULT_N 10000

using namespace MadVoro;

std::vector<Vector3D> GenerateRandomPoints(size_t N, const Vector3D &ll, const Vector3D &ur)
{
    std::vector<Vector3D> result;

    boost::random::uniform_real_distribution<double> distribution_x(ll.x, ur.x);
    boost::random::uniform_real_distribution<double> distribution_y(ll.y, ur.y);
    boost::random::uniform_real_distribution<double> distribution_z(ll.z, ur.z);

    boost::random::mt19937 rng;
    rng.seed(std::chrono::system_clock::now().time_since_epoch().count());

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

    size_t N = (argc == 2)? std::stoul(argv[1]) : DEFAULT_N;

    Vector3D ll(-1, -1, -1), ur(1, 1, 1);
    std::vector<Vector3D> points = GenerateRandomPoints(N, ll, ur);

    std::cout << "Generated " << N << " random points" << std::endl;
    std::cout << "Starting build" << std::endl;

    Voronoi3D voronoi(ll, ur);

    std::chrono::time_point<std::chrono::system_clock> start, end;

    start = std::chrono::system_clock::now();
    voronoi.Build(points);
    end = std::chrono::system_clock::now();


    double elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start).count();
    std::cout << "End of build, time taken is " << elapsed << " seconds" << std::endl;
    std::cout << "Rate (points/second, for a single processor): " << (points.size() / elapsed) << std::endl;

    return EXIT_SUCCESS;
}