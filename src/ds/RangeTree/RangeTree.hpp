/**
 * Implements a d-dimensional data structure for fast range queries, called range trees.
 * See more here: https://www.cs.umd.edu/class/fall2020/cmsc420-0201/Lects/lect17-range-tree.pdf
 * 
 * @author maor miz
*/
#ifndef RANGETREE_HPP
#define RANGETREE_HPP

#include "ds/BinaryTree.hpp"

namespace MadVoro
{
    namespace DataStructure
    {
        template<typename T>
        class RangeTree : public BinaryTree<T>
        {
        protected:
            using DefaultFilterFunction = std::function<bool(const T&)>;

            class RangeNode : public BinaryTree<T>::Node
            {
            public:
                RangeNode(): BinaryTree<T>::Node(), subtree(nullptr){};
                RangeNode(const T &value): BinaryTree<T>::Node(value), subtree(nullptr){};
                ~RangeNode() override{delete this->subtree;};

                RangeTree *subtree;
            };

            inline typename RangeTree<T>::RangeNode *getRoot() override{return this->actualRoot;};
            inline const typename RangeTree<T>::RangeNode *getRoot() const override{return this->actualRoot;};
            inline void setRoot(typename BinaryTree<T>::Node *other) override{if(other == nullptr) return; this->actualRoot = dynamic_cast<RangeNode*>(other); assert(this->actualRoot != nullptr);};

        public:
            RangeTree(RangeNode *root, const typename BinaryTree<T>::Compare &compare, int currentDim, int dimensions): BinaryTree<T>(root, compare), actualRoot(root), dim(dimensions), currentDim(currentDim){};
            RangeTree(const typename BinaryTree<T>::Compare &compare, int currentDim, int dimensions): RangeTree<T>(nullptr, compare, currentDim, dimensions){};
            RangeTree(int currentDim, int dimensions): RangeTree<T>([currentDim](const T &a, const T &b){return (a[currentDim] <= b[currentDim]);}, currentDim, dimensions){};
            RangeTree(int dimensions): RangeTree<T>(0, dimensions){};
            ~RangeTree() override{this->deleteSubtree(this->getRoot());};

            template<typename InputIterator>
            inline RangeTree(const InputIterator &first, const InputIterator &last, int dimensions): RangeTree<T>(dimensions){for(InputIterator it = first; it != last; it++) this->insert(*it);};

            template<typename RandomAccessIterator>
            inline void build(RandomAccessIterator first, RandomAccessIterator last){assert(this->treeSize == 0); std::sort(first, last, this->compare); this->setRoot(this->buildHelper(first, last));};
            template<typename Container>
            inline void build(Container &&container){this->build(container.begin(), container.end());};
            template<typename Container>
            inline void build(Container &container){this->build(container.begin(), container.end());};

            template<typename U = T, typename FilterFunction = DefaultFilterFunction>
            inline std::vector<T> range(const std::vector<std::pair<typename U::coord_type, typename U::coord_type>> &range, size_t maxToGet = std::numeric_limits<size_t>::max(), const FilterFunction &filter = [](const T&){return true;}) const
            {
                std::vector<T> result; 
                this->rangeHelper(range, this->getRoot(), 0, result, maxToGet, filter); 
                return result;
            };
            
            template<typename U = T, typename FilterFunction = DefaultFilterFunction>
            std::vector<T> circularRange(const U &center, typename U::coord_type, size_t maxToGet = std::numeric_limits<size_t>::max(), const FilterFunction &filter = [](const T&){return true;}) const;

        private:
            RangeNode *actualRoot = nullptr;
            int dim, currentDim;

            template<typename RandomAccessIterator>
            RangeNode *buildHelper(const RandomAccessIterator &first, const RandomAccessIterator &last);

            template<typename U = T, typename FilterFunction = DefaultFilterFunction>
            void rangeHelper(const std::vector<std::pair<typename U::coord_type, typename U::coord_type>> &range, const typename BinaryTree<T>::Node *node, int coord, std::vector<T> &result, size_t maxToGet, const FilterFunction &filter) const;
            template<typename U = T, typename FilterFunction = DefaultFilterFunction>
            void circularRangeHelper(const U &center, typename U::coord_type radius, const typename BinaryTree<T>::Node *node, int coord, std::vector<T> &result, size_t maxToGet, const FilterFunction &filter) const;
        };

        template<typename T>
        template<typename RandomAccessIterator>
        typename RangeTree<T>::RangeNode *RangeTree<T>::buildHelper(const RandomAccessIterator &first, const RandomAccessIterator &last)
        {
            if(first == last)
            {
                return nullptr;
            }
            RandomAccessIterator mid = first + (last - first) / 2;
            RangeNode *node = new RangeNode(*mid);
            ++this->treeSize;
            RangeNode *left = this->buildHelper(first, mid);
            node->left = left;
            RangeNode *right = this->buildHelper(mid + 1, last);
            node->right = right;
            
            if(left != nullptr) left->parent = node;
            if(right != nullptr) right->parent = node;

            this->updateNodeInfo(node);
            if(this->currentDim == this->dim - 1)
            {
                node->subtree = nullptr;
            }
            else
            {
                node->subtree = new RangeTree<T>(this->currentDim + 1, this->dim);
                node->subtree->build(this->getAllDecendants(node));
            }
            return node;
        }

        template<typename T>
        template<typename U, typename FilterFunction>
        void RangeTree<T>::rangeHelper(const std::vector<std::pair<typename U::coord_type, typename U::coord_type>> &range, const typename BinaryTree<T>::Node *node, int coord, std::vector<T> &result, size_t maxToGet, const FilterFunction &filter) const
        {
            if(node == nullptr or coord >= this->dim)
            {
                return;
            }
            if(node->minInSubtree[coord] > range[coord].second or node->maxInSubtree[coord] < range[coord].first)
            {
                return;
            }
            if(result.size() >= maxToGet)
            {
                return;
            }
            if(coord == this->dim - 1)
            {
                // search over the tree to find the maching points
                this->rangeHelper(range, node->left, coord, result, maxToGet, filter);
                if(node->value[coord] >= range[coord].first and node->value[coord] <= range[coord].second)
                {
                    if(filter(node->value))
                    {
                        result.push_back(node->value);
                    }
                }
                this->rangeHelper(range, node->right, coord, result, maxToGet, filter);
            }
            else
            {
                if(node->minInSubtree[coord] >= range[coord].first and node->maxInSubtree[coord] <= range[coord].second)
                {
                    // the whole subtree matches this coordinate! move on to the next one
                    assert(dynamic_cast<const RangeNode*>(node)->subtree != nullptr); // if the subtree was nullptr, then coord was dim-1.
                    dynamic_cast<const RangeNode*>(node)->subtree->rangeHelper(range, dynamic_cast<const RangeNode*>(node)->subtree->getRoot(), coord + 1, result, maxToGet, filter);
                }
                else
                {
                    // not the whole subtree matches, so check myself and move the recursion to my children
                    int i;
                    for(i = coord; i < this->dim; i++) if(node->value[i] < range[i].first or node->value[i] > range[i].second) break;
                    if(i == this->dim)
                    {
                        if(filter(node->value))
                        {
                            result.push_back(node->value);
                        }
                    }
                    this->rangeHelper(range, node->right, coord, result, maxToGet, filter);
                    this->rangeHelper(range, node->left, coord, result, maxToGet, filter);
                }
            }
        }

        template<typename T>
        template<typename U, typename FilterFunction>
        void RangeTree<T>::circularRangeHelper(const U &center, typename U::coord_type radius, const typename BinaryTree<T>::Node *node, int coord, std::vector<T> &result, size_t maxToGet, const FilterFunction &filter) const
        {
            if(node == nullptr or coord >= this->dim)
            {
                return;
            }
            
            if(result.size() >= maxToGet)
            {
                return;
            }

            typename U::coord_type cornerMax = center[coord] + radius;
            typename U::coord_type cornerMin = center[coord] - radius;

            if(node->minInSubtree[coord] > cornerMax or node->maxInSubtree[coord] < cornerMin)
            {
                return;
            }

            if(coord == this->dim - 1)
            {
                typename T::coord_type distanceSquared = typename T::coord_type();
                for(int i = 0; i < this->dim; i++)
                {
                    distanceSquared += ((node->value[i] - center[i]) * (node->value[i] - center[i]));
                }
                if(distanceSquared <= radius * radius)
                {
                    if(filter(node->value))
                    {
                        result.push_back(node->value);
                    }
                }
                if(node->minInSubtree[coord] <= cornerMax)
                {
                    this->circularRangeHelper(center, radius, node->left, coord, result, maxToGet, filter);
                }
                if(node->maxInSubtree[coord] >= cornerMin)
                {
                    this->circularRangeHelper(center, radius, node->right, coord, result, maxToGet, filter);
                }
            }
            else
            {
                if(node->minInSubtree[coord] >= cornerMin and node->maxInSubtree[coord] <= cornerMax)
                {
                    // the whole subtree matches this coordinate! move on to the next one
                    assert(dynamic_cast<const RangeNode*>(node)->subtree != nullptr); // if the subtree was nullptr, then coord was dim-1.
                    dynamic_cast<const RangeNode*>(node)->subtree->circularRangeHelper(center, radius, dynamic_cast<const RangeNode*>(node)->subtree->getRoot(), coord + 1, result, maxToGet, filter);
                }
                else
                {
                    typename T::coord_type distanceSquared = typename T::coord_type();
                    for(int i = 0; i < this->dim; i++)
                    {
                        distanceSquared += ((node->value[i] - center[i]) * (node->value[i] - center[i]));
                    }
                    if(distanceSquared <= (radius * radius))
                    {
                        if(filter(node->value))
                        {
                            result.push_back(node->value);
                        }
                    }
                    this->circularRangeHelper(center, radius, node->right, coord, result, maxToGet, filter);
                    this->circularRangeHelper(center, radius, node->left, coord, result, maxToGet, filter);
                }
            }
        }

        template<typename T>
        template<typename U, typename FilterFunction>
        inline std::vector<T> RangeTree<T>::circularRange(const U &center, typename U::coord_type radius, size_t maxToGet, const FilterFunction &filter) const
        {
            // a simple implementation - 
            std::vector<T> result;
            std::vector<std::pair<typename U::coord_type, typename U::coord_type>> rectangle;
            for(int i = 0; i < this->dim; i++)
            {
                rectangle.push_back(std::make_pair<typename U::coord_type, typename U::coord_type>(center[i] - radius, center[i] + radius));
            }

            for(const auto &point : this->range(rectangle, maxToGet, filter))
            {
                typename T::coord_type squaredDistance = typename T::coord_type();
                for(int i = 0; i < this->dim; i++)
                {
                    squaredDistance += (center[i] - point[i]) * (center[i] - point[i]);
                }
                if(squaredDistance <= radius * radius)
                {
                    result.push_back(point);
                }
            }
            return result;
        }
    }
}

#endif // RANGETREE_HPP