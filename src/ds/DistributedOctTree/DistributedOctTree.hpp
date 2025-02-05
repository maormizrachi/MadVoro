#ifndef DISTRIBUTED_OCTTREE_HPP
#define DISTRIBUTED_OCTTREE_HPP

#ifdef MADVORO_WITH_MPI

#include <vector>
#include <assert.h>
#include <utility>
#include <array>
#include <bitset>
#include <boost/container/flat_set.hpp>
#include <mpi.h>
#include "ds/OctTree/OctTree.hpp"
#include "RankedValue.hpp"

#ifdef DEBUG_MODE
#include <iostream>
#endif // DEBUG_MODE

#define MAX_DIRECTIONS_SIZE (MAX_DEPTH + 1)

#define MAX_RANKS_FOR_LEAF_DEFAULT 4

namespace MadVoro
{    
    namespace DataStructure
    {
        template<typename T>
        using BoundingBox = MadVoro::Geometry::BoundingBox<T>;
        
        template<typename T>
        using Sphere = MadVoro::Geometry::Sphere<T>;

        template<typename T, int max_ranks_per_leaf, int max_directions_size>
        using DistributedOctNode = typename OctTree<OctTreeRankedValue<T, max_ranks_per_leaf, max_directions_size>>::OctTreeNode;

        template<typename T, int max_ranks_per_leaf = MAX_RANKS_FOR_LEAF_DEFAULT>
        class DistributedOctTree
        {
        public:
            using RankedValue = OctTreeRankedValue<T, max_ranks_per_leaf, MAX_DIRECTIONS_SIZE>;
            using DistributedOctTreeNode = DistributedOctNode<T, max_ranks_per_leaf, MAX_DIRECTIONS_SIZE>;

        public:
            DistributedOctTree(const OctTree<T> *tree, bool detailedNodeInfo = false, const MPI_Comm &comm = MPI_COMM_WORLD);

            ~DistributedOctTree(){delete this->octTree;};

            #ifdef DEBUG_MODE
            void print() const{this->octTree->print();};
            #endif // DEBUG_MODE

            boost::container::flat_set<int> getIntersectingRanks(const Sphere<T> &sphere) const;
            
            inline boost::container::flat_set<int> getIntersectingRanks(const T &center, const typename T::coord_type radius) const{return this->getIntersectingRanks(Sphere<T>(center, radius));};
            
            int getDepth() const{return this->octTree->getDepth();};
            
            OctTree<RankedValue> *getOctTree(){return this->octTree;};

            const OctTree<RankedValue> *getOctTree() const{return this->octTree;};

            std::vector<T> getRankValues(int _rank) const;

            inline std::vector<T> getMyValues() const{return this->getRankValues(this->rank);};

            std::vector<BoundingBox<T>> getRankBoundingBoxes(int _rank) const;

            std::vector<BoundingBox<T>> getRankBoundingBoxesUnique(int _rank) const;

            inline std::vector<BoundingBox<T>> getMyBoundingBoxes() const{return this->getRankBoundingBoxes(this->rank);};

            inline std::vector<BoundingBox<T>> getMyBoundingBoxesUnique() const{return this->getRankBoundingBoxesUnique(this->rank);};

            std::vector<std::vector<BoundingBox<T>>> getBoundingBoxesOfRanks(void) const;

            std::vector<BoundingBox<T>> getBigBoundingBoxesOfRanks(void) const;

            std::vector<std::vector<direction_t>> getRankDirections(int _rank) const;

            inline std::vector<const DistributedOctTreeNode*> getRankNodes(int _rank) const
            {
                return this->getValuesIf([](const DistributedOctTreeNode *node){return node->value.owners.empty();}, 
                                    [&_rank](const DistributedOctTreeNode *node){return std::find(node->value.owners.cbegin(), node->value.owners.cend(), _rank) != node->value.owners.cend();});
            }

            inline std::vector<std::vector<direction_t>> getMyDirections() const{return this->getRankDirections(this->rank);};

            inline std::vector<const DistributedOctTreeNode*> getMyNodes() const{return this->getRankNodes(this->rank);};

            template<typename U>
            std::vector<std::pair<typename T::coord_type, typename T::coord_type>> getClosestFurthestPointsByRanks(const U &point) const;

            #ifdef DEBUG_MODE
            inline bool validate() const{if(this->octTree != nullptr) return this->validateHelper(this->octTree->getRoot()); return true;};
            #endif // DEBUG_MODE

            std::vector<const DistributedOctTreeNode*> getValuesIf(const std::function<bool(const DistributedOctTreeNode*)> ifOpenFunction, const std::function<bool(const DistributedOctTreeNode*)> &ifAddValueFunction) const;

        private:
            OctTree<RankedValue> *octTree = nullptr;
            MPI_Comm comm;
            int rank, size;
            bool detailedNodeInfo; // whether or not to save detailed info on the leaf nodes
            size_t treeSize;

            void buildTreeHelper(DistributedOctTreeNode *newNode, const typename OctTree<T>::OctTreeNode *node, std::vector<direction_t> &directionsInMyTree);

            void buildTree(const OctTree<T> *tree);

            #ifdef DEBUG_MODE
            bool validateHelper(const DistributedOctTreeNode *node) const;
            #endif // DEBUG_MODE
        };

        template<typename T, int max_ranks_per_leaf>
        DistributedOctTree<T, max_ranks_per_leaf>::DistributedOctTree(const OctTree<T> *tree, bool detailedNodeInfo, const MPI_Comm &comm): comm(comm), detailedNodeInfo(detailedNodeInfo)
        {
            MPI_Comm_rank(this->comm, &this->rank);
            MPI_Comm_size(this->comm, &this->size);
            this->buildTree(tree);
        }

        template<typename T, int max_ranks_per_leaf>
        void DistributedOctTree<T, max_ranks_per_leaf>::buildTreeHelper(DistributedOctTreeNode *newNode, const typename OctTree<T>::OctTreeNode *node, std::vector<direction_t> &directionsInMyTree)
        {
            assert(newNode != nullptr);
            unsigned char valueToSend = 0; // assumes `CHILDREN` is 8. this variable contains 1 in the `i`th bit iff child `i` exists
            if(node != nullptr)
            {
                for(int i = 0; i < CHILDREN; i++)
                {
                    bool bit = (node->children[i] != nullptr || (node->isLeaf and newNode->getChildNumberContaining(node->value) == i));
                    valueToSend |= (bit << i);
                }
            }

            std::vector<unsigned char> childBuff(this->size);
            MPI_Allgather(&valueToSend, 1, MPI_UNSIGNED_CHAR, &childBuff[0], 1, MPI_BYTE, this->comm);

            for(int i = 0; i < CHILDREN; i++)
            {
                bool recursiveBuild = false;
                int containingValue = UNDEFINED_OWNER; // who has the `i`th child
                std::vector<int> ranks_containing_child;

                for(int _rank = 0; _rank < this->size; _rank++)
                {
                    if(((childBuff[_rank] >> i) & 0x1) == 1)
                    {
                        ranks_containing_child.push_back(_rank);

                        if(containingValue == UNDEFINED_OWNER)
                        {
                            // first rank to announce it holds the value
                            containingValue = _rank;
                        }
                        else
                        {
                            // more than one child has a point in this route of the tree, so we continue to split recursively
                            recursiveBuild = true;
                        }
                    }
                }
                // if somebody holds the value
                if(containingValue != UNDEFINED_OWNER)
                {
                    // someone holds the `i`th child
                    this->treeSize++;
                    newNode->createChild(i); // creates the child in my own tree
                    if(recursiveBuild and (ranks_containing_child.size() > max_ranks_per_leaf))
                    {
                        // there are several holders, call recursive build (until we reach one holder)
                        // determine what's the next node in my own tree to continue the recursive build with
                        // this node might be null, if I don't have any nodes this depth in the tree
                        const typename OctTree<T>::OctTreeNode *nextNode = nullptr;

                        bool advancedNode = false; // if went down to a child of current node
                        if(node != nullptr)
                        {
                            if(node->isLeaf)
                            {
                                nextNode = newNode->children[i]->boundingBox.contains(node->value)? node : nullptr;
                            }
                            else
                            {
                                advancedNode = true;
                                directionsInMyTree.push_back(i);
                                nextNode = node->children[i];
                            }
                        }
                        else
                        {
                            nextNode = nullptr;
                        }
                        // continue recursively
                        this->buildTreeHelper(newNode->children[i], nextNode, directionsInMyTree);
                        if(advancedNode)
                        {
                            directionsInMyTree.pop_back();
                        }
                        // several owners
                        newNode->children[i]->value.owners.resize(0);
                        newNode->children[i]->isLeaf = false;
                    }
                    else
                    {
                        if((ranks_containing_child.size() == 1) or (this->detailedNodeInfo))
                        {
                            // there is only one holder, set its owner field
                            if(this->detailedNodeInfo)
                            {
                                if(this->rank == containingValue)
                                {
                                    const typename OctTree<T>::OctTreeNode *childNode;
                                    if(node->isLeaf)
                                    {
                                        childNode = node;
                                    }
                                    else
                                    {
                                        childNode = node->children[i];
                                        directionsInMyTree.push_back(i);
                                    }
                    
                                    // copy value, and directions
                                    std::memcpy(&newNode->children[i]->value.value, &childNode->value, sizeof(T));
                                    std::memcpy(newNode->children[i]->value.directions, directionsInMyTree.data(), sizeof(direction_t) * directionsInMyTree.size());
                                    newNode->children[i]->value.directions[directionsInMyTree.size()] = PATH_END_DIRECTION;

                                    if(!node->isLeaf)
                                    {
                                        directionsInMyTree.pop_back();
                                    }
                                }
                                MPI_Bcast(&newNode->children[i]->value.value, sizeof(T), MPI_BYTE, containingValue, this->comm);
                                MPI_Bcast(newNode->children[i]->value.directions, sizeof(direction_t) * MAX_DIRECTIONS_SIZE, MPI_BYTE, containingValue, this->comm);
                            }

                            newNode->children[i]->value.owners = {containingValue};
                        }
                        else
                        {
                            // several owners
                            newNode->children[i]->value.owners.resize(ranks_containing_child.size());
                            std::copy(ranks_containing_child.cbegin(), ranks_containing_child.cend(), newNode->children[i]->value.owners.begin());
                        }
                        newNode->children[i]->isLeaf = true;
                    }
                }
                else
                {
                    this->treeSize++;
                    newNode->children[i] = nullptr;
                }
            }
        }

        #ifdef DEBUG_MODE
        template<typename T, int max_ranks_per_leaf>
        bool DistributedOctTree<T, max_ranks_per_leaf>::validateHelper(const DistributedOctTreeNode *node) const
        {
            if(node == nullptr)
            {
                return true;
            }
            bool hasChildren = false;
            for(int i = 0; i < CHILDREN; i++)
            {
                if(node->children[i] != nullptr)
                {
                    hasChildren = true;
                    break;
                }
            }
            if(hasChildren)
            {
                assert(node->isLeaf);
            }
            for(int i = 0; i < CHILDREN; i++)
            {
                assert(this->validateHelper(node->children[i]));
            }
            return true;
        }
        #endif // DEBUG_MODE

        template<typename T, int max_ranks_per_leaf>
        void DistributedOctTree<T, max_ranks_per_leaf>::buildTree(const OctTree<T> *tree)
        {
            assert(this->octTree == nullptr);
            if(tree == nullptr or tree->getRoot() == nullptr)
            {
                return;
            }
            std::vector<direction_t> directions;
            directions.reserve(MAX_DIRECTIONS_SIZE);
            // tree->print();

            RankedValue ll, ur;

            this->octTree = new OctTree<RankedValue>(tree->getLL(), tree->getUR());
            this->treeSize = 0;
            this->buildTreeHelper(this->octTree->getRoot(), tree->getRoot(), directions);
        }

        template<typename T, int max_ranks_per_leaf>
        boost::container::flat_set<int> DistributedOctTree<T, max_ranks_per_leaf>::getIntersectingRanks(const Sphere<T> &sphere) const
        {
            boost::container::flat_set<int> ranks;
            for(const RankedValue &point : this->octTree->range(Sphere<RankedValue>(sphere.center, sphere.radius)))
            {
                size_t numOwners = point.owners.size();
                for(size_t i = 0; i < numOwners; i++)
                {
                    ranks.insert(point.owners[i]);
                }
            }
            return ranks;
        }

        template<typename T, int max_ranks_per_leaf>
        template<typename U>
        std::vector<std::pair<typename T::coord_type, typename T::coord_type>> DistributedOctTree<T, max_ranks_per_leaf>::getClosestFurthestPointsByRanks(const U &point) const
        {
            const typename T::coord_type &maxVal = std::numeric_limits<typename T::coord_type>::max();
            const typename T::coord_type &minVal = std::numeric_limits<typename T::coord_type>::lowest();
            
            std::pair<T, T> initialPair = std::make_pair<T, T>(T(maxVal, maxVal, maxVal), T(minVal, minVal, minVal));
            std::vector<std::pair<typename T::coord_type, typename T::coord_type>> distances(this->size, {maxVal, minVal});

            std::vector<const DistributedOctTreeNode*> nodes;
            nodes.reserve(this->getDepth() * CHILDREN);
            nodes.push_back(this->octTree->getRoot());

            T closestPoint, furthestPoint;
            while(!nodes.empty())
            {
                const DistributedOctTreeNode *node = nodes.back();
                nodes.pop_back();

                if(node == nullptr)
                {
                    continue;
                }
                if(!node->isLeaf)
                {
                    for(int i = 0; i < CHILDREN; i++)
                    {
                        nodes.push_back(node->children[i]);
                    }
                    continue;
                }
                // node is a value node
                closestPoint = node->boundingBox.closestPoint(point);
                furthestPoint = node->boundingBox.furthestPoint(point);
                typename T::coord_type closestDist = 0, furthestDist = 0;
                #ifdef MADVORO_WITH_VCL
                    Vec8d diff(closestPoint[0] - point[0], closestPoint[1] - point[1], closestPoint[2] - point[2],
                                furthestPoint[0] - point[0], furthestPoint[1] - point[1], furthestPoint[2] - point[2],
                                0, 0);
                    Vec8d dist = diff * diff;
                    closestDist = dist[0] + dist[1] + dist[2];
                    furthestDist = dist[3] + dist[4] + dist[5];
                #else // MADVORO_WITH_VCL
                    closestDist += (closestPoint[0] - point[0]) * (closestPoint[0] - point[0]);
                    furthestDist += (furthestPoint[0] - point[0]) * (furthestPoint[0] - point[0]);
                    closestDist += (closestPoint[1] - point[1]) * (closestPoint[1] - point[1]);
                    furthestDist += (furthestPoint[1] - point[1]) * (furthestPoint[1] - point[1]);
                    closestDist += (closestPoint[2] - point[2]) * (closestPoint[2] - point[2]);
                    furthestDist += (furthestPoint[2] - point[2]) * (furthestPoint[2] - point[2]);
                #endif // MADVORO_WITH_VCL
                size_t numOwners = node->value.owners.size();

                for(size_t i = 0; i < numOwners; i++)
                {
                    int owner = node->value.owners[i];
                    if(distances[owner].first > closestDist)
                    {
                        distances[owner].first = closestDist;
                    }
                    if(distances[owner].second < furthestDist)
                    {
                        distances[owner].second = furthestDist;
                    }
                }
            }
            return distances;
        }


        template<typename T, int max_ranks_per_leaf>
        std::vector<T> DistributedOctTree<T, max_ranks_per_leaf>::getRankValues(int _rank) const
        {
            std::vector<T> values;
            std::vector<const DistributedOctTreeNode*> nodes;
            nodes.reserve(this->getDepth() * CHILDREN);

            nodes.push_back(this->octTree->getRoot());

            while(!nodes.empty())
            {
                const DistributedOctTreeNode *node = nodes.back();
                nodes.pop_back();

                if(node == nullptr)
                {
                    continue;
                }

                size_t numOwners = node->value.owners.size();

                if(numOwners == 0)
                {
                    for(int i = 0; i < CHILDREN; i++)
                    {
                        nodes.push_back(node->children[i]);
                    }
                }
                else
                {
                    if(std::find(node->value.owners.cbegin(), node->value.owners.cend(), _rank) != node->value.owners.cend())
                    {
                        // is an owner
                        values.emplace_back(node->value.value);
                    }
                }
            }
            return values;
        }

        template<typename T, int max_ranks_per_leaf>
        std::vector<BoundingBox<T>> DistributedOctTree<T, max_ranks_per_leaf>::getRankBoundingBoxes(int _rank) const
        {
            std::vector<BoundingBox<T>> boxes;
            std::vector<const DistributedOctTreeNode*> nodes;
            nodes.reserve(this->getDepth() * CHILDREN);

            nodes.push_back(this->octTree->getRoot());

            while(!nodes.empty())
            {
                const DistributedOctTreeNode *node = nodes.back();
                nodes.pop_back();

                if(node == nullptr)
                {
                    continue;
                }

                size_t numOwners = node->value.owners.size();

                if(numOwners == 0)
                {
                    // several owners
                    for(int i = 0; i < CHILDREN; i++)
                    {
                        nodes.push_back(node->children[i]);
                    }
                }
                else
                {
                    if(std::find(node->value.owners.cbegin(), node->value.owners.cend(), _rank) != node->value.owners.cend())
                    {
                        // is an owner
                        boxes.emplace_back(BoundingBox<T>(node->boundingBox.getLL(), node->boundingBox.getUR()));
                    }
                }
            }
            return boxes;
        }

        template<typename T, int max_ranks_per_leaf>
        std::vector<BoundingBox<T>> DistributedOctTree<T, max_ranks_per_leaf>::getRankBoundingBoxesUnique(int _rank) const
        {
            std::vector<BoundingBox<T>> boxes;
            std::vector<const DistributedOctTreeNode*> nodes;
            nodes.reserve(this->getDepth() * CHILDREN);

            nodes.push_back(this->octTree->getRoot());

            while(!nodes.empty())
            {
                const DistributedOctTreeNode *node = nodes.back();
                nodes.pop_back();

                if(node == nullptr)
                {
                    continue;
                }

                size_t numOwners = node->value.owners.size();

                if(numOwners == 0)
                {
                    for(int i = 0; i < CHILDREN; i++)
                    {
                        nodes.push_back(node->children[i]);
                    }
                }
                else
                {
                    if(*std::min_element(node->value.owners.cbegin(), node->value.owners.cend()) == _rank)
                    {
                        // is an owner
                        boxes.emplace_back(BoundingBox<T>(node->boundingBox.getLL().value, node->boundingBox.getUR().value));
                    }
                }
            }
            return boxes;
        }

        template<typename T, int max_ranks_per_leaf>
        std::vector<std::vector<BoundingBox<T>>> DistributedOctTree<T, max_ranks_per_leaf>::getBoundingBoxesOfRanks(void) const
        {
            std::vector<std::vector<BoundingBox<T>>> boxes(this->size);
            std::vector<const DistributedOctTreeNode*> nodes;
            nodes.reserve(this->getDepth() * CHILDREN);

            nodes.push_back(this->octTree->getRoot());

            while(!nodes.empty())
            {
                const DistributedOctTreeNode *node = nodes.back();
                nodes.pop_back();

                if(node == nullptr)
                {
                    continue;
                }

                size_t numOwners = node->value.owners.size();

                if(numOwners == 0)
                {
                    for(int i = 0; i < CHILDREN; i++)
                    {
                        nodes.push_back(node->children[i]);
                    }
                }
                else
                {
                    for(int _rank : node->value.owners)
                    {
                        boxes[_rank].push_back(node->boundingBox);
                    }
                }
            }
            return boxes;
        }

        template<typename T, int max_ranks_per_leaf>
        std::vector<BoundingBox<T>> DistributedOctTree<T, max_ranks_per_leaf>::getBigBoundingBoxesOfRanks(void) const
        {
            const BoundingBox<T> &rootBB = this->octTree->getRoot()->boundingBox;
            std::vector<std::pair<T, T>> boxes(this->size, {rootBB.getUR(), rootBB.getLL()});
            std::vector<const DistributedOctTreeNode*> nodes;
            nodes.reserve(this->getDepth() * CHILDREN);

            nodes.push_back(this->octTree->getRoot());

            while(!nodes.empty())
            {
                const DistributedOctTreeNode *node = nodes.back();
                nodes.pop_back();

                if(node == nullptr)
                {
                    continue;
                }

                size_t numOwners = node->value.owners.size();

                if(numOwners == 0)
                {
                    for(int i = 0; i < CHILDREN; i++)
                    {
                        nodes.push_back(node->children[i]);
                    }
                }
                else
                {
                    const T &BB_ll = node->boundingBox.getLL();
                    const T &BB_ur = node->boundingBox.getUR();
                    for(int _rank : node->value.owners)
                    {
                        boxes[_rank].first[0] = std::min(boxes[_rank].first[0], BB_ll[0]);
                        boxes[_rank].second[0] = std::max(boxes[_rank].second[0], BB_ur[0]);
                        boxes[_rank].first[1] = std::min(boxes[_rank].first[1], BB_ll[1]);
                        boxes[_rank].second[1] = std::max(boxes[_rank].second[1], BB_ur[1]);
                        boxes[_rank].first[2] = std::min(boxes[_rank].first[2], BB_ll[2]);
                        boxes[_rank].second[2] = std::max(boxes[_rank].second[2], BB_ur[2]);
                    }
                }
            }
            std::vector<BoundingBox<T>> returnBoundingBoxes;
            for(int _rank = 0; _rank < this->size; _rank++)
            {
                returnBoundingBoxes.push_back(BoundingBox<T>(boxes[_rank].first, boxes[_rank].second));
            }
            return returnBoundingBoxes;
        }


        template<typename T, int max_ranks_per_leaf>
        std::vector<std::vector<direction_t>> DistributedOctTree<T, max_ranks_per_leaf>::getRankDirections(int _rank) const
        {
            std::vector<std::vector<direction_t>> directions;
            std::vector<const DistributedOctTreeNode*> nodes;
            nodes.reserve(this->getDepth() * CHILDREN);

            nodes.push_back(this->octTree->getRoot());

            while(!nodes.empty())
            {
                const DistributedOctTreeNode *node = nodes.back();
                nodes.pop_back();

                if(node == nullptr)
                {
                    continue;
                }

                size_t numOwners = node->value.owners.size();

                if(numOwners == 0)
                {
                    for(int i = 0; i < CHILDREN; i++)
                    {
                        nodes.push_back(node->children[i]);
                    }
                }
                else
                {
                    if(std::find(node->value.owners.cbegin(), node->value.owners.cend(), _rank) != node->value.owners.cend())
                    {
                        // is an owner
                        size_t directionsLength = std::distance(node->value.directions, std::find(node->value.directions, node->value.directions + MAX_DIRECTIONS_SIZE, PATH_END_DIRECTION)) + 1; // + 1 for the `PATH_END_DIRECTION`
                        directions.push_back(std::vector(node->value.directions, node->value.directions + directionsLength));
                    }
                }
            }
            return directions;
        }

        template<typename T, int max_ranks_per_leaf>
        std::vector<const typename DistributedOctTree<T, max_ranks_per_leaf>::DistributedOctTreeNode*>
            DistributedOctTree<T, max_ranks_per_leaf>::getValuesIf(const std::function<bool(const DistributedOctTreeNode*)> ifOpenFunction, const std::function<bool(const DistributedOctTreeNode*)> &ifAddValueFunction) const
        {
            std::vector<const DistributedOctTreeNode*> nodes = {this->octTree->getRoot()};
            nodes.reserve(this->getDepth() * CHILDREN);

            std::vector<const DistributedOctTreeNode*> result;

            while(not nodes.empty())
            {
                const DistributedOctTreeNode *node = nodes.back();
                nodes.pop_back();
                if(node == nullptr)
                {
                    continue;
                }
                if(node->isLeaf)
                {
                    if(ifAddValueFunction(node))
                    {
                        result.push_back(node);
                    }
                    continue;
                }
                if(ifOpenFunction(node))
                {
                    for(int i = 0; i < CHILDREN; i++)
                    {
                        nodes.push_back(node->children[i]);
                    }
                }
            }
            return result;
        }
    }
}

#endif // MADVORO_WITH_MPI

#endif // DISTRIBUTED_OCTTREE_HPP