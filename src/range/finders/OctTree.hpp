#ifndef _OCT_TREE_FINDER_HPP
#define _OCT_TREE_FINDER_HPP

#include "ds/OctTree/OctTree.hpp"
#include "utils/IndexedVector.hpp"
#include "RangeFinder.hpp"

#define DIMENSIONS 3

namespace MadVoro
{
    namespace Range
    {
        class OctTreeFinder : public RangeFinder
        {
        public:
            template<typename RandomAccessIterator>
            OctTreeFinder(RandomAccessIterator first, RandomAccessIterator last, const Point3D &ll ,const Point3D &ur);
            
            inline OctTreeFinder(std::vector<Point3D> &myPoints, const Point3D &ll ,const Point3D &ur): OctTreeFinder(myPoints.begin(), myPoints.end(), ll, ur){};
            
            inline OctTreeFinder(MadVoro::DataStructure::OctTree<IndexedPoint3D> *tree, const std::vector<Point3D> &myPoints): octTree(tree), myPoints(myPoints), givenOctTree(true){}

            inline ~OctTreeFinder() override
            {
                if(not this->givenOctTree)
                {
                    delete this->octTree;
                }
            };

            std::vector<size_t> closestPointInSphere(const Point3D &center, double radius, const Point3D &point, const _set<size_t> &ignore) const override
            {
                std::pair<IndexedPoint3D, double> closestPointPair = this->octTree->getClosestPointInSphere(Sphere<Point3D>(center, radius), point,
                                                                                                            [&ignore](const IndexedPoint3D &vec){return ignore.find(vec.getIndex()) == ignore.cend();});
                const IndexedPoint3D &closestPoint = closestPointPair.first;
                const double &closestDistance = closestPointPair.second;

                if(closestDistance != std::numeric_limits<typename IndexedPoint3D::coord_type>::max())
                {
                    return std::vector<size_t>({closestPoint.index});
                }
                return std::vector<size_t>(); // empty
            };

            inline std::vector<size_t> range(const Point3D &center, double radius, size_t N, const _set<size_t> &ignore) const override
            {
                std::vector<size_t> toReturn;
                for(const IndexedPoint3D &vec : this->octTree->range(Sphere<IndexedPoint3D>(IndexedPoint3D(center.x, center.y, center.z, ILLEGAL_IDX), radius + EPSILON), N,
                                                                    [&ignore](const IndexedPoint3D &vec){return ignore.find(vec.getIndex()) == ignore.cend();}))
                {
                    toReturn.push_back(vec.index);
                }
                return toReturn;
            };

            inline size_t size() const override{return this->octTree->getSize();};

            inline const Point3D &getPoint(size_t index) const override{return this->myPoints[index];};

        private:

            std::vector<Point3D> myPoints;
            MadVoro::DataStructure::OctTree<IndexedPoint3D> *octTree;
            bool givenOctTree;
        };

        template<typename RandomAccessIterator>
        inline OctTreeFinder::OctTreeFinder(RandomAccessIterator first, RandomAccessIterator last, const Point3D &ll ,const Point3D &ur)
        {
            this->octTree = new MadVoro::DataStructure::OctTree<IndexedPoint3D>(ll, ur);
            this->givenOctTree = false;
            this->myPoints.reserve(last - first);
            size_t index = 0;
            for(RandomAccessIterator it = first; it != last; it++)
            {
                const Point3D &vec = *it;
                this->myPoints.push_back(vec);
                IndexedPoint3D idx_vec = IndexedPoint3D(vec.x, vec.y, vec.z, index);
                this->octTree->insert(idx_vec);
                index++;
            }
        }
    }
}

#endif // _OCT_TREE_FINDER_HPP