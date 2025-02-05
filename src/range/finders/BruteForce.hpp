#ifndef _BRUTE_FORCE_RANGE_HPP
#define _BRUTE_FORCE_RANGE_HPP

#include "RangeFinder.hpp"

namespace MadVoro
{
    namespace Range
    {
        class BruteForceFinder : public RangeFinder
        {
        public:
            template<typename RandomAccessIterator>
            BruteForceFinder(const RandomAccessIterator &first, const RandomAccessIterator &last):
                        points(std::vector<Point3D>()){this->points.insert(this->points.begin(), first, last); this->pointsSize = this->points.size();};
            
            BruteForceFinder(const std::vector<Point3D> &points): BruteForceFinder(points.begin(), points.end()){};

            inline ~BruteForceFinder() = default;
            
            std::vector<size_t> closestPointInSphere(const Point3D &center, double radius, const Point3D &point, const _set<size_t> &ignore) const override
            {
                size_t closestSoFarIndex = std::numeric_limits<size_t>::max();
                double closestSoFar = std::numeric_limits<double>::max();
                const Point3D *_points = this->points.data();
                for(size_t i = 0; i < this->pointsSize; i++)
                {
                    if((ignore.find(i) != ignore.cend()))
                    {
                        continue;
                    }
                    //  __builtin_prefetch(&_points[i]);
                    const Point3D &_point = _points[i];
                    double sphere_dx = _point.x - center.x;
                    double sphere_dy = _point.y - center.y;
                    double sphere_dz = _point.z - center.z;
                    double dx = point.x - center.x;
                    double dy = point.y - center.y;
                    double dz = point.z - center.z;
                    sphere_dx *= sphere_dx;
                    sphere_dy *= sphere_dy;
                    sphere_dz *= sphere_dz;
                    dx *= dx;
                    dy *= dy;
                    dz *= dz;
                    double distanceToPoint = dx + dy + dz;
                    double distanceToSphere = sphere_dx + sphere_dy + sphere_dz;
                    if((distanceToPoint < closestSoFar) and (distanceToSphere <= ((radius * radius) + EPSILON)))
                    {
                        closestSoFar = distanceToPoint;
                        closestSoFarIndex = i;
                    }
                }

                std::vector<size_t> result;
                result.reserve(1);
                if(closestSoFarIndex != std::numeric_limits<size_t>::max())
                {
                    // a point was found
                    result.push_back(closestSoFarIndex);
                }
                return result;
            }

            inline const Point3D &getPoint(size_t index) const override{return this->points[index];};

            inline std::vector<size_t> range(const Point3D &center, double radius, size_t N, const _set<size_t> &ignore) const override
            {
                std::vector<size_t> result;
                const Point3D *_points = this->points.data();
                for(size_t i = 0; i < this->pointsSize; i++)
                {
                    if(result.size() >= N)
                    {
                        break;
                    }

                    if((ignore.find(i) != ignore.cend()))
                    {
                        continue;
                    }

                    const Point3D &point = _points[i];
                    double dx = point.x - center.x;
                    double dy = point.y - center.y;
                    double dz = point.z - center.z;
                    dx *= dx;
                    dy *= dy;
                    dz *= dz;
                    if((dx + dy + dz) <= ((radius * radius) + EPSILON))
                    {
                        result.push_back(i);
                    }
                }
                return result;
            }

            inline size_t size() const override{return this->pointsSize;};

        private:
            std::vector<Point3D> points;
            size_t pointsSize;
        };
    }
}

#endif // _BRUTE_FORCE_RANGE_HPP
