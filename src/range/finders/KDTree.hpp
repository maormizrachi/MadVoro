#ifndef _KDTREE_FINDER_HPP
#define _KDTREE_FINDER_HPP

#include "ds/KDTree/KDTree.hpp"
#include "utils/IndexedVector.hpp"
#include "RangeFinder.hpp"

#define DIMENSIONS 3

namespace MadVoro
{
    namespace Range
    {
        class KDTreeFinder : public RangeFinder
        {
        public:
            template<typename RandomAccessIterator>
            KDTreeFinder(RandomAccessIterator first, RandomAccessIterator last, const Point3D &ll ,const Point3D &ur);
            
            inline KDTreeFinder(std::vector<Point3D> &myPoints, const Point3D &ll ,const Point3D &ur): KDTreeFinder(myPoints.begin(), myPoints.end(), ll, ur){};
            
            inline ~KDTreeFinder() override{delete this->kdTree;};
            
            std::vector<size_t> closestPointInSphere(const Point3D &center, double radius, const Point3D &point, const _set<size_t> &ignore) const override
            {
                std::pair<IndexedPoint3D, double> closestPointPair = this->kdTree->getClosestPointInSphere(Sphere<Point3D>(center, radius), point,
                                                                                                            [&ignore](const IndexedPoint3D &vec){return ignore.find(vec.getIndex()) == ignore.cend();});
                const IndexedPoint3D &closestPoint = closestPointPair.first;
                const double &closestDistance = closestPointPair.second;

                if(closestDistance != std::numeric_limits<typename IndexedPoint3D::coord_type>::max())
                {
                    return std::vector<size_t>({closestPoint.index});
                }
                return std::vector<size_t>(); // empty
            }

            inline const Point3D &getPoint(size_t index) const override{return this->myPoints[index];};

            inline std::vector<size_t> range(const Point3D &center, double radius, size_t N, const _set<size_t> &ignore) const override
            {
                std::vector<size_t> toReturn;
                for(const IndexedPoint3D &vec : this->kdTree->range(Sphere<Point3D>(center, radius), N,
                                                                    [&ignore](const IndexedPoint3D &vec){return ignore.find(vec.getIndex()) == ignore.cend();}))
                {
                    toReturn.push_back(vec.index);
                }
                return toReturn;
            };
            inline size_t size() const override{return this->kdTree->getSize();};

        private:
            std::vector<Point3D> myPoints;
            MadVoro::DataStructure::KDTree<IndexedPoint3D, DIMENSIONS> *kdTree;
        };

        template<typename RandomAccessIterator>
        inline KDTreeFinder::KDTreeFinder(RandomAccessIterator first, RandomAccessIterator last, const Point3D &ll ,const Point3D &ur)
        {
            size_t index = 0;
            this->kdTree = new MadVoro::DataStructure::KDTree<IndexedPoint3D, DIMENSIONS>(ll, ur);

            myPoints.reserve(last - first);
            for(RandomAccessIterator it = first; it != last; it++)
            {
                const Point3D &vec = *it;
                IndexedPoint3D idx_vec = IndexedPoint3D(vec.x, vec.y, vec.z, index);
                this->myPoints.push_back(vec);
                this->kdTree->insert(idx_vec);
                index++;
            }
        }
    }
}

#endif // _KDTREE_FINDER_HPP