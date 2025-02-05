#ifndef _RANGE_TREE_FINDER_HPP
#define _RANGE_TREE_FINDER_HPP

#include "ds/BinaryTree.hpp"
#include "ds/RangeTree/RangeTree.hpp"
#include "utils/IndexedVector.hpp"
#include "RangeFinder.hpp"

#define DIMENSIONS 3

namespace MadVoro
{
    namespace Range
    {
        class RangeTreeFinder : public RangeFinder
        {
        public:
            template<typename RandomAccessIterator>
            RangeTreeFinder(RandomAccessIterator first, RandomAccessIterator last);
            
            inline RangeTreeFinder(std::vector<Point3D> &myPoints): RangeTreeFinder(myPoints.begin(), myPoints.end()){};
            
            inline ~RangeTreeFinder() override{delete this->rangeTree;};

            inline const Point3D &getPoint(size_t index) const override{return this->myPoints[index];};

            inline std::vector<size_t> range(const Point3D &center, double radius, size_t N, const _set<size_t> &ignore) const override
            {
                std::vector<size_t> toReturn;
                for(const IndexedPoint3D &vec : this->rangeTree->circularRange(center, radius, N, [&ignore](const IndexedPoint3D &vec){return ignore.find(vec.getIndex()) == ignore.cend();}))
                {
                    toReturn.push_back(vec.index);
                }
                return toReturn;
            };
            inline size_t size() const override{return this->rangeTree->size();};

        private:
            std::vector<Point3D> myPoints;
            MadVoro::DataStructure::RangeTree<MadVoro::Range::IndexedPoint3D> *rangeTree;
        };

        template<typename RandomAccessIterator>
        inline MadVoro::Range::RangeTreeFinder::RangeTreeFinder(RandomAccessIterator first, RandomAccessIterator last)
        {
            std::vector<IndexedPoint3D> data;
            size_t index = 0;

            myPoints.reserve(last - first);
            for(RandomAccessIterator it = first; it != last; it++)
            {
                const Point3D &vec = *it;
                this->myPoints.push_back(vec);
                IndexedPoint3D idx_vec = IndexedPoint3D(vec.x, vec.y, vec.z, index);
                data.push_back(idx_vec);
                index++;
            }
            this->rangeTree = new MadVoro::DataStructure::RangeTree<IndexedPoint3D>(DIMENSIONS);
            this->rangeTree->build(data.begin(), data.end());
        }
    }
}

#endif // _RANGE_TREE_FINDER_HPP