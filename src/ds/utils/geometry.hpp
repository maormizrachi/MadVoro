#ifndef GEOMETRY_UTILS_HPP
#define GEOMETRY_UTILS_HPP

#include "BoundingBox.hpp"
#include "Sphere.hpp"

#define TOLERANCE 1e-12

namespace MadVoro
{
    namespace Geometry
    {
        /**
         * @brief Checks if a sphere intersects a box
        */
        template<typename T, typename U>
        bool SphereBoxIntersection(const BoundingBox<T> &box, const Sphere<U> &sphere);
    }
}

template<typename T, typename U>
bool MadVoro::Geometry::SphereBoxIntersection(const BoundingBox<T> &box, const Sphere<U> &sphere)
{
    typename T::coord_type distance = 0;
    const T &boxLL = box.getLL(), &boxUR = box.getUR();
    for(int i = 0; i < DIM; i++)
    {
        typename T::coord_type centerCoord = sphere.center[i];
        typename T::coord_type closestPointCoord;

        if(centerCoord < boxLL[i])
        {
            closestPointCoord = boxLL[i];
        }
        else
        {
            if(centerCoord <= boxUR[i])
            {
                closestPointCoord = centerCoord;
            }
            else
            {
                closestPointCoord = boxUR[i];
            }
        }
        typename T::coord_type _distance = (closestPointCoord - sphere.center[i]);
        _distance *= _distance;
        distance += _distance;
    }
    return (distance <= ((1 + TOLERANCE) * (sphere.radius * sphere.radius)));
};

#endif // GEOMETRY_UTILS_HPP