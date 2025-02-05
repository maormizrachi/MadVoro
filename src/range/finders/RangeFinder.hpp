#ifndef _RANGE_FINDER_HPP
#define _RANGE_FINDER_HPP

#include <boost/container/flat_set.hpp>
#include <boost/unordered_set.hpp>
#include <vector>
#include <limits>
#include <unordered_set>
#include "elementary/Point3D.hpp"

namespace MadVoro
{
    namespace Range
    {
        template<typename T>
        using Sphere = MadVoro::Geometry::Sphere<T>;

        template<typename T>
        using _set = boost::container::flat_set<T>; // boost::unordered_set<T>; // std::unordered_set<T>; // boost::container::flat_set<T>;
    }
}

namespace MadVoro
{
    namespace Range
    {
        class RangeFinder
        {
        public:
            virtual ~RangeFinder() = default;
            
            virtual std::vector<size_t> range(const Point3D &center, double radius, size_t N = std::numeric_limits<size_t>::max(), const _set<size_t> &ignore = _set<size_t>()) const
            {
                throw MadVoro::Exception::MadVoroException("RangeFinder::range: method not implemented");
            }
            
            virtual std::vector<size_t> closestPointInSphere(const Point3D &center, double radius, const Point3D &point, const _set<size_t> &ignore) const
            {
                throw MadVoro::Exception::MadVoroException("RangeFinder::closestPointInSphere: method not implemented");
            }
            
            virtual const Point3D &getPoint(size_t index) const = 0;
            
            virtual size_t size() const = 0;
        };
    }
}

#endif // _RANGE_FINDER_HPP