#ifndef GROUP_RANGETREE_HPP
#define GROUP_RANGETREE_HPP

#include <vector>
#include <algorithm>
#include "exception/MadVoroException.hpp"
#include "ds/GroupTree/GroupTree.hpp"

namespace MadVoro
{
    namespace DataStructure
    {
        template<typename T, int N>
        class GroupRangeTree : public GroupTree<T, N>
        {
        protected:
            using DefaultFilterFunction = std::function<bool(const T&)>;

            class GroupRangeNode : public GroupTree<T, N>::Node
            {
            public:
                GroupRangeNode(): GroupTree<T, N>::Node(), subtree(nullptr){};
                GroupRangeNode(const T &value): GroupTree<T, N>::Node(value), subtree(nullptr){};
                template<typename Container>
                GroupRangeNode(const Container &container): GroupTree<T, N>::Node(container), subtree(nullptr){};
                ~GroupRangeNode() override{delete this->subtree;};

                GroupRangeTree<T, N> *subtree;
            };

            void splitNode(typename GroupTree<T, N>::Node *node) override;
            GroupRangeNode *tryInsert(const T &value);

        private:
            inline typename GroupRangeTree<T, N>::GroupRangeNode *getRoot() override{return this->actualRoot;};
            inline const typename GroupRangeTree<T, N>::GroupRangeNode *getRoot() const override{return this->actualRoot;};
            inline void setRoot(typename GroupTree<T, N>::Node *other) override{if(other == nullptr) return; this->actualRoot = dynamic_cast<GroupRangeNode*>(other); assert(this->actualRoot != nullptr);};

            template<typename RandomAccessIterator>
            GroupRangeNode *fastBuildHelper(const RandomAccessIterator &first, const RandomAccessIterator &last);
            template<typename U, typename FilterFunction>
            void rangeHelper(const std::vector<std::pair<typename U::coord_type, typename U::coord_type>> &range, const typename GroupTree<T, N>::Node *node, int coord, std::vector<T> &result, size_t maxToGet, const FilterFunction &filter) const;

            void recreateSubtree(GroupRangeNode *node);

            inline void deleteHelper(typename GroupTree<T, N>::Node *node) override{if(node == nullptr) return; deleteHelper(node->left); deleteHelper(node->right); delete dynamic_cast<GroupRangeNode*>(node);};
            inline void deleteTree() override{this->deleteHelper(this->getRoot()); this->setRoot(nullptr); this->treeSize = 0;};

        public:
            GroupRangeTree(GroupRangeNode *root, const typename GroupTree<T, N>::Compare &compare, int currentDim, int dimensions): GroupTree<T, N>(root, compare), actualRoot(root), dim(dimensions), currentDim(currentDim){};
            GroupRangeTree(const typename GroupTree<T, N>::Compare &compare, int currentDim, int dimensions): GroupRangeTree(nullptr, compare, currentDim, dimensions){};
            GroupRangeTree(int currentDim, int dimensions): GroupRangeTree([currentDim](const T &a, const T &b){return (a[currentDim] <= b[currentDim]);}, currentDim, dimensions){};
            GroupRangeTree(int dimensions): GroupRangeTree<T, N>(0, dimensions){};
            template<typename InputIterator>
            inline GroupRangeTree(const InputIterator &first, const InputIterator &last, int dimensions): GroupRangeNode(){for(InputIterator it = first; it != last; it++) this->insert(*it);};
            ~GroupRangeTree() override{this->deleteTree();};

            template<typename U = T, typename FilterFunction = DefaultFilterFunction>
            inline std::vector<T> range(const std::vector<std::pair<typename U::coord_type, typename U::coord_type>> &range, size_t maxToGet = std::numeric_limits<size_t>::max(), const FilterFunction &filter = [](const T&){return true;}) const
            {
                std::vector<T> result;
                this->rangeHelper(range, this->getRoot(), 0, result, maxToGet, filter);
                return result;
            };
            template<typename U = T, typename FilterFunction = DefaultFilterFunction>
            std::vector<T> circularRange(const U &center, typename U::coord_type radius, size_t maxToGet = std::numeric_limits<size_t>::max(), const FilterFunction &filter = [](const T&){return true;}) const; 
            
            template<typename RandomAccessIterator>
            inline void build(RandomAccessIterator first, RandomAccessIterator last);

            template<typename Container>
            inline void build(Container &&container){this->build(container.begin(), container.end());};
            
            template<typename Container>
            inline void build(Container &container){this->build(container.begin(), container.end());};

        private:
            GroupRangeNode *actualRoot = nullptr;
            int dim, currentDim;
        };
    }
}

template<typename T, int N>
void MadVoro::DataStructure::GroupRangeTree<T, N>::recreateSubtree(GroupRangeNode *node)
{
    if(node == nullptr)
    {
        return;
    }
    if(this->currentDim == this->dim - 1)
    {
        return;
    }
    delete node->subtree;
    node->subtree = new GroupRangeTree<T, N>(this->currentDim + 1, this->dim);
    node->subtree->build(this->getAllDecendants(node));
    
    // rebuild also my children:
    this->recreateSubtree(dynamic_cast<GroupRangeNode*>(node->left));
    this->recreateSubtree(dynamic_cast<GroupRangeNode*>(node->right));
}

template<typename T, int N>
template<typename RandomAccessIterator>
typename MadVoro::DataStructure::GroupRangeTree<T, N>::GroupRangeNode *MadVoro::DataStructure::GroupRangeTree<T, N>::fastBuildHelper(const RandomAccessIterator &first, const RandomAccessIterator &last)
{
    if(first == last)
    {
        return nullptr;
    }

    GroupRangeNode *node;
    if(last - first <= N)
    {
        node = new GroupRangeNode();
        node->numValues = last - first;
        int j = 0;
        for(RandomAccessIterator it = first; it != last; it++)
        {
            node->values[j++] = *it;
        }
    }
    else
    {
        RandomAccessIterator mid = first + (last - first) / 2;
        node = new GroupRangeNode(*mid);
        GroupRangeNode *left = this->fastBuildHelper(first, mid);
        node->left = left;
        GroupRangeNode *right = this->fastBuildHelper(mid + 1, last); 
        node->right = right;
    }
    if(this->currentDim == this->dim - 1)
    {
        node->subtree = nullptr;
    }
    else
    {
        node->subtree = new GroupRangeTree<T, N>(this->currentDim + 1, this->dim);
        node->subtree->build(this->getAllDecendants(node));
    }
    this->treeSize += node->numValues;
    
    if(node->left != nullptr)
    {
        node->left->parent = node;
        this->updateNodeInfo(node->left);
    }
    if(node->right != nullptr)
    {
        node->right->parent = node;
        this->updateNodeInfo(node->right);
    }
    return node;
}

template<typename T, int N>
template<typename RandomAccessIterator>
inline void MadVoro::DataStructure::GroupRangeTree<T, N>::build(RandomAccessIterator first, RandomAccessIterator last)
{
    if(this->treeSize != 0)
    {
        throw MadVoro::Exception::MadVoroException("Can not build a GroupRangeTree, because tree is not empty");
    }
    std::sort(first, last, this->compare);
    this->setRoot(this->fastBuildHelper(first, last));
    this->updateNodeInfo(this->getRoot());
};

template<typename T, int N>
void MadVoro::DataStructure::GroupRangeTree<T, N>::splitNode(typename GroupTree<T, N>::Node *node)
{
    if(node == nullptr)
    {
        return;
    }
    // assert(node->numValues > 1);
    
    std::sort(node->values.begin(), node->values.begin() + node->numValues, this->compare);
    int numLeft = node->numValues / 2;
    int numRight = node->numValues - numLeft - 1;
    std::array<T, N> left;
    std::array<T, N> right;
    for(int i = 0; i < numLeft; i++)
    {
        left[i] = node->values[i];
    }
    for(int i = numLeft + 1; i < node->numValues; i++)
    {
        right[i - (numLeft + 1)] = node->values[i];
    }

    if(numLeft > 0 and node->left == nullptr)
    {
        node->left = new GroupRangeNode(left);
        node->left->numValues = numLeft;
        node->left->parent = node;
    }
    if(numRight > 0 and node->right == nullptr)
    {
        node->right = new GroupRangeNode(right);
        node->right->numValues = numRight;
        node->right->parent = node;
    }

    node->numValues = 1;
    node->values[0] = node->values[numLeft];
    node->minInSubtree = node->maxInSubtree = node->values[0];

    this->updateNodeInfo(node->left);
    this->updateNodeInfo(node->right);
    this->recreateSubtree(dynamic_cast<GroupRangeNode*>(node)); // rebuild the sub-trees
}

template<typename T, int N>
typename MadVoro::DataStructure::GroupRangeTree<T, N>::GroupRangeNode *MadVoro::DataStructure::GroupRangeTree<T, N>::tryInsert(const T &value)
{
    if(this->getRoot() == nullptr)
    {
        this->setRoot(new GroupRangeNode(value));
        if(this->currentDim != this->dim - 1)
        {
            this->getRoot()->subtree = new GroupRangeTree<T, N>(this->currentDim + 1, this->dim);
            this->getRoot()->subtree->insert(value);
        }
        return this->getRoot();
    }
    GroupRangeNode *current = this->getRoot();
    while(true)
    {
        if(this->compare(value, current->values[0]))
        {
            if(current->left == nullptr)
            {
                break;
            }
            current = dynamic_cast<GroupRangeNode*>(current->left);
        }
        else
        {
            if(current->right == nullptr)
            {
                break;
            }
            current = dynamic_cast<GroupRangeNode*>(current->right);
        }
    }

    // check if value already exists:
    for(int i = 0; i < current->numValues; i++)
    {
        if(value == current->values[i])
        {
            return current;
        }
    }

    GroupRangeNode *addTo;
    if(current->numValues < N and current->left == nullptr and current->right == nullptr)
    {
        addTo = current;
    }
    else
    {
        if(current->numValues < N)
        {
            this->splitNode(current);
        }
        if(this->compare(value, current->values[0]))
        {
            if(current->left == nullptr)
            {
                // not a leaf, so we need to create it a child
                current->left = new GroupRangeNode(value);
                current->left->parent = current;
                current->left->numValues = 0;
            }
            addTo = dynamic_cast<GroupRangeNode*>(current->left);
        }
        else
        {
            if(current->right == nullptr)
            {
                // not a leaf, so we need to create it a child
                current->right = new GroupRangeNode(value);
                current->right->parent = current;
                current->right->numValues = 0;
            }
            addTo = dynamic_cast<GroupRangeNode*>(current->right);
        }
    }
    addTo->values[addTo->numValues++] = value;

    // add the new value to its ancestors subtrees:
    if(this->currentDim != this->dim - 1)
    {
        GroupRangeNode *node = addTo;
        while(node != nullptr)
        {
            if(node->subtree == nullptr)
            {
                assert(node == addTo);
                // should not happen, unless `node` is `addTo` (because its ancesors' subtrees should exist)
                node->subtree = new GroupRangeTree<T, N>(this->currentDim + 1, this->dim);
            }
            node->subtree->insert(value);
            node = dynamic_cast<GroupRangeNode*>(node->parent);
        }
    }
    return addTo;
}

template<typename T, int N>
template<typename U, typename FilterFunction>
void MadVoro::DataStructure::GroupRangeTree<T, N>::rangeHelper(const std::vector<std::pair<typename U::coord_type, typename U::coord_type>> &range, const typename GroupTree<T, N>::Node *node, int coord, std::vector<T> &result, size_t maxToGet, const FilterFunction &filter) const
{
    if(node == nullptr or coord >= this->dim)
    {
        return;
    }
    if((node->minInSubtree[coord] > range[coord].second) or (node->maxInSubtree[coord] < range[coord].first))
    {
        return;
    }
    if(result.size() >= maxToGet)
    {
        return; // don't continue
    }
    if(coord == this->dim - 1)
    {
        // search over the tree to find the maching points
        this->rangeHelper(range, node->left, coord, result, maxToGet, filter);
        for(int i = 0; i < node->numValues; i++)
        {
            const T &value = node->values[i];
            if((value[coord] >= range[coord].first) and (value[coord] <= range[coord].second) and (result.size() < maxToGet))
            {
                if(filter(value))
                {
                    result.push_back(value);
                }
            }
        }
        this->rangeHelper(range, node->right, coord, result, maxToGet, filter);
    }
    else
    {
        if(node->minInSubtree[coord] >= range[coord].first and node->maxInSubtree[coord] <= range[coord].second)
        {
            // the whole subtree matches this coordinate! move on to the next one
            assert(dynamic_cast<const GroupRangeNode*>(node)->subtree != nullptr); // if the subtree was nullptr, then coord was dim-1.
            dynamic_cast<const GroupRangeNode*>(node)->subtree->rangeHelper(range, dynamic_cast<const GroupRangeNode*>(node)->subtree->getRoot(), coord + 1, result);
        }
        else
        {
            for(int i = 0; i < node->numValues; i++)
            {
                const T &value = node->values[i];
                int j;
                for(j = coord; j < this->dim; j++) if(value[j] < range[j].first or value[j] > range[j].second) break;
                if((j == this->dim) and (result.size() < maxToGet))
                {
                    if(filter(value))
                    {
                        result.push_back(value);
                    }
                }
            }
            this->rangeHelper(range, node->right, coord, result, maxToGet, filter);
            this->rangeHelper(range, node->left, coord, result, maxToGet, filter);
        }
    }
}

template<typename T, int N>
template<typename U, typename FilterFunction>
std::vector<T> MadVoro::DataStructure::GroupRangeTree<T, N>::circularRange(const U &center, typename U::coord_type radius, size_t maxToGet, const FilterFunction &filter) const
{
    std::vector<std::pair<typename U::coord_type, typename U::coord_type>> rectangularRange;
    for(int i = 0; i < this->dim; i++)
    {
        rectangularRange.push_back(std::make_pair(center[i] - radius, center[i] + radius));
    }
    std::vector<T> result;
    for(const T &possible_value : this->range(rectangularRange, maxToGet))
    {
        typename T::coord_type distanceSquared = 0;
        for(int i = 0; i < this->dim; i++)
        {
            distanceSquared += (possible_value[i] - center[i]) * (possible_value[i] - center[i]);
        }
        if(distanceSquared <= (radius * radius))
        {
            result.push_back(possible_value);
        }
    }
    return result;
}

#endif // GROUP_RANGETREE_HPP