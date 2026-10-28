// The box cells are built in: its faces, walls and periodic axes
#ifndef MADVORO_BUILDER_BUILD_DOMAIN_HPP
#define MADVORO_BUILDER_BUILD_DOMAIN_HPP

#include <array>
#include <cstddef>
#include <vector>

#include "elementary/Face3D.hpp"
#include "utils/Predicates3D.hpp"

namespace MadVoro
{

template <typename PointT>
bool PointInPoly(const std::vector<Face3D<PointT>> &faces, const PointT &point)
{
    std::array<PointT, 4> vec;
    vec[3] = point;
    for(const Face3D<PointT> &face : faces)
    {
        vec[0] = face.vertices[0];
        vec[1] = face.vertices[1];
        vec[2] = face.vertices[2];
        if(orient3d(vec) > 0)
        {
            return false;
        }
    }
    return true;
}

template <typename PointT>
bool PointInDomain(const PointT &ll, const PointT &ur, const PointT &point)
{
    return point.x > ll.x and point.x < ur.x and point.y > ll.y and point.y < ur.y and point.z > ll.z and point.z < ur.z;
}

//! \brief The translation of a periodic image from the box [ll, ur) (zero inside it)
template <typename PointT>
PointT PeriodicTranslationOf(const PointT &p, const PointT &ll, const PointT &ur, const std::array<bool, 3> &periodic)
{
    PointT translation;
    const PointT L = ur - ll;
    if(periodic[0])
    {
        if(p.x < ll.x)
        {
            translation.x = -L.x;
        }
        else if(p.x >= ur.x)
        {
            translation.x = L.x;
        }
    }
    if(periodic[1])
    {
        if(p.y < ll.y)
        {
            translation.y = -L.y;
        }
        else if(p.y >= ur.y)
        {
            translation.y = L.y;
        }
    }
    if(periodic[2])
    {
        if(p.z < ll.z)
        {
            translation.z = -L.z;
        }
        else if(p.z >= ur.z)
        {
            translation.z = L.z;
        }
    }
    return translation;
}

template <typename PointT>
struct BuildDomain
{
    PointT ll, ur;
    std::vector<Face3D<PointT>> boxFaces; // empty: the box [ll, ur]
    std::vector<Face3D<PointT>> box; // the faces mirrors are made across
    std::vector<PointT> normals; // of `box`
    std::array<bool, 3> periodic = {false, false, false};

    bool HasPeriodicBoundaries(void) const
    {
        return this->periodic[0] or this->periodic[1] or this->periodic[2];
    }

    bool IsOutside(const PointT &point) const
    {
        return this->boxFaces.empty() ? not PointInDomain(this->ll, this->ur, point) : not PointInPoly(this->boxFaces, point);
    }

    PointT PeriodicTranslation(const PointT &p) const
    {
        return PeriodicTranslationOf(p, this->ll, this->ur, this->periodic);
    }
};

} // namespace MadVoro

#endif // MADVORO_BUILDER_BUILD_DOMAIN_HPP
