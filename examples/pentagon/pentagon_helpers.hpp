#ifndef PENTAGON_HELPERS
#define PENTAGON_HELPERS

#include <vector>
#include <random>
#include <algorithm>
#include <mpi.h>
#include "../Vector3D.hpp"
#include <Voronoi3D.hpp>
#include <elementary/PointOps.hpp>
using namespace MadVoro::fallback;

using Hyperplane3D = std::pair<Vector3D, double>;

std::vector<Hyperplane3D> GetHyperplanes(const std::vector<MadVoro::Face<Vector3D>> &faces, const Vector3D &pointInPoly)
{
    std::vector<Hyperplane3D> hyperplanes;
    for(const auto &face : faces)
    {
        Vector3D normal = normalize(CrossProduct(face.vertices[1] - face.vertices[0], face.vertices[2] - face.vertices[1]));
        double D = ScalarProd(normal, face.vertices[0]);
        if(ScalarProd(normal, pointInPoly) > D)
        {
            normal *= -1;
            D *= -1;
        }
        hyperplanes.emplace_back(normal, D);
    }
    return hyperplanes;
}

bool PointInPolygon(const std::vector<Hyperplane3D> &planes, const Vector3D &point)
{
    for(const Hyperplane3D &plane : planes)
    {
        if(ScalarProd(plane.first, point) > plane.second)
        {
            return false;
        }
    }
    return true;
}

std::vector<Vector3D> GeneratePoints(const std::vector<MadVoro::Face<Vector3D>> &faces, size_t N, const Vector3D &ll, const Vector3D &ur)
{
    std::vector<std::uniform_real_distribution<double>> distributions;
    for(int i = 0; i < 3; i++)
    {
        distributions.push_back(std::uniform_real_distribution(ll[i], ur[i]));
    }
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    std::mt19937 re(rank);
    std::vector<Vector3D> points;
    Vector3D pointInPoly = (faces[0].vertices[0] + faces[1].vertices[2]) / 2;
    std::vector<Hyperplane3D> planes = GetHyperplanes(faces, pointInPoly); 

    points.reserve(N);
    while(points.size() != N)
    {
        Vector3D point;
        for(int i = 0; i < 3; i++)
        {
            point[i] = distributions[i](re);
        }

        if(PointInPolygon(planes, point) and std::find(points.begin(), points.end(), point) == points.end())
        {
            points.push_back(point);
        }
    }
    return points;
}

#endif // PENTAGON_HELPERS
