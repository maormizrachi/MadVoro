#ifndef _OCTTREE_HPP
#define _OCTTREE_HPP

#include <limits>
#include <vector>
#include <functional>
#include <assert.h>
#include <utility>
#include <type_traits>
#include <stack>

#ifdef DEBUG_MODE
#include <iostream>
#endif // DEBUG_MODE

#include "ds/utils/raw_type.h"
#include "ds/utils/geometry.hpp"
#include "exception/MadVoroException.hpp"

#define DIM 3
#define CHILDREN 8 // 2^DIM
#define PATH_END_DIRECTION (-1)
#define MAX_DEPTH 64

namespace MadVoro
{
    namespace DataStructure
    {
        template<typename T>
        using BoundingBox = MadVoro::Geometry::BoundingBox<T>;

        template<typename T>
        using Sphere = MadVoro::Geometry::Sphere<T>;

        typedef size_t octnode_id_t;
        typedef char direction_t;

        template<typename T>
        class OctTree
        {
            template<typename U>
            class DistributedOctTree;
            template<typename U>
            friend class DistributedOctTree;

            using Raw_type = typename is_raw_type_defined<T>::type;
            using DefaultFilterFunction = std::function<bool(const T&)>;

        public:
            class OctTreeNode
            {
                friend class OctTree;

            public:
                inline OctTreeNode(const T &ll, const T &ur): isLeaf(false), value((ll + ur)/2), boundingBox(BoundingBox<Raw_type>(ll, ur)), parent(nullptr), height(0), depth(0)
                {
                    for(int i = 0; i < CHILDREN; i++)
                    {
                        this->children[i] = nullptr;
                    }
                }

                inline OctTreeNode(const T &point): isLeaf(true), value(point), boundingBox(BoundingBox<Raw_type>(point, point)), parent(nullptr), height(0), depth(0)
                {
                    for(int i = 0; i < CHILDREN; i++)
                    {
                        this->children[i] = nullptr;
                    }
                }

                inline OctTreeNode(OctTreeNode &&other): isLeaf(other.isLeaf), value(other.value), boundingBox(other.boundingBox), height(other.height), depth(other.depth)
                {
                    for(int i = 0; i < CHILDREN; i++)
                    {
                        this->children[i] = other.children[i];
                        other.children[i] = nullptr;
                    }
                    this->parent = other.parent;
                    other.parent = nullptr;
                }

                OctTreeNode(OctTreeNode *parent, int childNumber);

                virtual ~OctTreeNode() = default;
                
                template<typename U>
                int getChildNumberContaining(const U &point) const;

                template<typename U>
                const OctTreeNode *getChildContaining(const U &point) const{return this->children[this->getChildNumberContaining(point)];};

                virtual OctTreeNode *addLeafChild(int childIndex, const T &point);

                virtual OctTreeNode *createChild(int childNumber);

                #ifdef DEBUG_MODE
                inline void print() const
                {
                    std::cout << this->value << ", BB: " << this->boundingBox.getLL() << ", " << this->boundingBox.getUR() << " (depth: " << this->depth << ", height: " << this->height << ")" << std::endl;
                }
                #endif // DEBUG_MODE

                bool isLeaf; // if a leaf
                T value; // if a leaf, that's a point value, otherwise, thats the value for partition
                BoundingBox<Raw_type> boundingBox; // the bounding box this node induces
                std::array<OctTreeNode*, CHILDREN> children; // if a leaf, all children are nullptr
                OctTreeNode *parent;
                int height; // height of a leaf is 0
                int depth; // depth of the root is 0
                octnode_id_t id;

            protected:
                void fixHeightsRecursively();
                void splitNode();
            };

        protected:
            void deleteSubtree(OctTreeNode *node);

            template<typename U>
            const OctTreeNode *tryFind(const U &point) const;
            
            template<typename U>
            inline OctTreeNode *tryFind(const U &point){return const_cast<OctTreeNode*>(std::as_const(*this).tryFind(point));};

            template<typename U>
            const OctTreeNode *tryFindParent(const U &point) const;

            template<typename U>
            inline OctTreeNode *tryFindParent(const U &point){return const_cast<OctTreeNode*>(std::as_const(*this).tryFindParent(point));};

            template<typename U>
            OctTreeNode *tryInsert(const U &point);

            #ifdef DEBUG_MODE
            void printHelper(const OctTreeNode *node, int indent) const;
            #endif // DEBUG_MODE

            void getAllDecendantsHelper(const OctTreeNode *node, std::vector<T> &result) const;

            std::vector<T> getAllDecendants(const OctTreeNode *node) const;
            
            OctTreeNode *root;
            T ll, ur;
            size_t treeSize;
            size_t nodesNumber;
            mutable std::vector<const OctTreeNode*> nodes_stack;

        public:
            explicit OctTree(const T &ll, const T &ur): root(nullptr), treeSize(0), nodesNumber(0){this->setBounds(ll, ur);};

            template<typename InputIterator>
            explicit OctTree(const T &ll, const T &ur, const InputIterator &first, const InputIterator &last);

            template<typename Container>
            inline OctTree(const T &ll, const T &ur, Container container): OctTree(ll, ur, container.begin(), container.end()){};

            inline explicit OctTree(): root(nullptr), treeSize(0), nodesNumber(0){};

            virtual inline ~OctTree(){this->deleteSubtree(this->getRoot());};

            inline void clear()
            {
                this->deleteSubtree(this->getRoot());
                this->setRoot(nullptr);
                this->treeSize = 0;
                this->nodesNumber = 0;
            }
            
            template<typename U>
            inline bool insert(const U &point)
            {
                OctTreeNode *newNode = this->tryInsert(point);
                if(newNode != nullptr)
                {
                    newNode->id = (this->nodesNumber++);
                    this->treeSize++;
                    return true;
                }
                return false;
            };

            template<typename U>
            inline bool find(const U &point) const{return this->tryFind(point) != nullptr;};

            template<typename U>
            T &findValue(const U &point);

            template<typename U>
            OctTreeNode *findNode(const U &point);

            template<typename U>
            const OctTreeNode *findNodeContainingBoundingBox(const BoundingBox<U> &boundingBox) const; // TODO: necessary? What about just find value of that node?

            template<typename U>
            inline T findParent(const U &point) const
            {
                const OctTreeNode *parent = this->tryFindParent(point);
                if(parent == nullptr)
                {
                    MadVoro::Exception::MadVoroException eo("OctTree: could not find a parent (is the tree empty or point out of boundaries?)");
                    eo.addEntry("Value", point);
                    throw eo;
                }
                return parent->value;
            };

            template<typename U>
            inline bool remove(const U &point)
            {
                OctTreeNode *node = this->tryFind(point);
                if(node == nullptr)
                {
                    return false;
                }
                OctTreeNode *newParent = this->removeLeaf(node);
                if(newParent != nullptr)
                {
                    newParent->fixHeightsRecursively();
                }
                this->treeSize--;
                return true;
            }

            OctTreeNode *removeLeaf(OctTreeNode *node);

            virtual inline OctTreeNode *getRoot(){return this->root;};

            virtual inline const OctTreeNode *getRoot() const{return this->root;};

            virtual inline void setRoot(OctTreeNode *other){this->root = other;};

            virtual void setBounds(const T &ll, const T &ur)
            {
                this->ll = ll;
                this->ur = ur;
                assert(this->getRoot() == nullptr);
                this->setRoot(new OctTreeNode(ll, ur));
                this->getRoot()->parent = nullptr;
            }

            const T &getLL() const{return this->ll;};
            const T &getUR() const{return this->ur;};

            #ifdef DEBUG_MODE
            void print() const{this->printHelper(this->getRoot(), 0);};
            #endif // DEBUG_MODE

            inline int getDepth() const
            {
                if(this->getRoot() == nullptr)
                {
                    throw MadVoro::Exception::MadVoroException("OctTree::getDepth: root is null");
                }
                return this->getRoot()->height;
            };

            inline size_t getSize() const{return this->treeSize;};

            template<typename U, typename FilterFunction = DefaultFilterFunction>
            std::vector<T> range(const Sphere<U> &sphere, size_t N = std::numeric_limits<size_t>::max(), const FilterFunction &filter = [](const T&){return true;}) const;

            template<typename U, typename FilterFunction = DefaultFilterFunction>
            std::pair<T, typename T::coord_type> getClosestPointInSphere(const Sphere<U> &sphere, const T &point, const FilterFunction &filter = [](const T&){return true;}) const;

            const OctTreeNode *getNodeByDirections(const direction_t *directions = nullptr) const;

            template<typename U>
            std::pair<T, typename T::coord_type> getClosestPointInfo(const U &point, bool includeSelf = true) const;

            template<typename U>
            inline T closestPoint(const U &point, bool includeSelf = true) const{return this->getClosestPointInfo(point, includeSelf).first;};

            template<typename U>
            inline typename T::coord_type closestPointDistance(const U &point, bool includeSelf = true) const{return this->getClosestPointInfo(point, includeSelf).second;};
        };

        template<typename T>
        void OctTree<T>::deleteSubtree(OctTreeNode *node)
        {
            if(node == nullptr)
            {
                return;
            }
            for(int i = 0; i < CHILDREN; i++)
            {
                this->deleteSubtree(node->children[i]);
            }
            delete node;
        }

        template<typename T>
        OctTree<T>::OctTreeNode::OctTreeNode(OctTreeNode *parent, int childNumber): isLeaf(false), parent(parent), height(0), depth(0)
        {
            assert(parent != nullptr);

            // determine box:
            Raw_type new_ll, new_ur;
            const Raw_type &parentLL = parent->boundingBox.getLL(), &parentUR = parent->boundingBox.getUR();
            for(int i = 0; i < DIM; i++)
            {
                if((childNumber >> ((DIM - 1) - i)) & 1)
                {
                    new_ll[i] = (parentLL[i] + parentUR[i]) / 2;
                    new_ur[i] = parentUR[i];
                }
                else
                {
                    new_ll[i] = parentLL[i];
                    new_ur[i] = (parentLL[i] + parentUR[i]) / 2;
                }
                this->value[i] = (new_ll[i] + new_ur[i]) / 2;
            }

            this->boundingBox.setBounds(new_ll, new_ur);

            for(int i = 0; i < CHILDREN; i++)
            {
                this->children[i] = nullptr;
            }
            this->fixHeightsRecursively();
        }

        template<typename T>
        typename OctTree<T>::OctTreeNode *OctTree<T>::OctTreeNode::addLeafChild(int childIndex, const T &point)
        {
            this->children[childIndex] = new OctTreeNode(point);
            this->children[childIndex]->parent = this;
            this->children[childIndex]->fixHeightsRecursively();
            return this->children[childIndex];
        }

        template<typename T>
        typename OctTree<T>::OctTreeNode *OctTree<T>::OctTreeNode::createChild(int childNumber)
        {
            assert(this->children[childNumber] == nullptr);
            this->children[childNumber] = new OctTreeNode(this, childNumber);
            return this->children[childNumber];
        }

        template<typename T>
        template<typename U>
        int OctTree<T>::OctTreeNode::getChildNumberContaining(const U &point) const
        {
            assert(this->boundingBox.contains(point));
            int direction = 0;
            for(int i = 0; i < DIM; i++)
            {
                direction = (direction << 1) | ((this->value[i] < point[i])? 1 : 0);
            }
            return direction;
        }

        template<typename T>
        template<typename U>
        const typename OctTree<T>::OctTreeNode *OctTree<T>::tryFindParent(const U &point) const
        {
            const OctTreeNode *current = this->getRoot();
            while(current != nullptr)
            {
                if(current->isLeaf)
                {
                    return current;
                }
                // otherwise, determine the direction to go
                const OctTreeNode *nextChild = current->getChildContaining(point);
                if(nextChild == nullptr)
                {
                    return current;
                }
                current = nextChild;
            }
            return nullptr;
        }

        template<typename T>
        template<typename U>
        const typename OctTree<T>::OctTreeNode *OctTree<T>::tryFind(const U &point) const
        {
            const OctTreeNode *current = this->getRoot();
            while(current != nullptr)
            {
                if(current->isLeaf)
                {
                    if(current->value == point)
                    {
                        return current;
                    }
                    return nullptr;
                }
                // otherwise, determine the direction to go
                current = current->getChildContaining(point);
            }
            return nullptr;
        }

        template<typename T>
        template<typename U>
        const typename OctTree<T>::OctTreeNode *OctTree<T>::findNodeContainingBoundingBox(const BoundingBox<U> &boundingBox) const
        {
            const OctTreeNode *current = this->getRoot();
            const U &BB_center = (boundingBox.getLL() + boundingBox.getUR()) * 0.5;
            while(current != nullptr)
            {
                if(boundingBox.contained(current->boundingBox))
                {
                    return current;
                }
                // otherwise, determine the direction to go
                current = current->getChildContaining(BB_center);
            }
            return nullptr;
        }

        template<typename T>
        void OctTree<T>::OctTreeNode::fixHeightsRecursively()
        {
            if(this->isLeaf)
            {
                this->height = 0;
            }
            if(this->parent != nullptr)
            {
                this->parent->height = std::max<int>(this->parent->height, this->height + 1);
                this->parent->fixHeightsRecursively();
            }
            this->depth = (this->parent == nullptr)? 0 : this->parent->depth + 1;
        }

        template<typename T>
        void OctTree<T>::OctTreeNode::splitNode()
        {
            assert(this->parent != nullptr);
            assert(this->isLeaf);

            // get my index
            int i;
            for(i = 0; i < CHILDREN; i++)
            {
                if(this == this->parent->children[i])
                {
                    break;
                }
            }

            // this node is the `i`th child of its parent
            // replace it with a new (non-value) node, which will be our parent
            this->parent->children[i] = nullptr;
            this->parent->createChild(i);
            
            this->parent = this->parent->children[i];
            this->parent->isLeaf = false;
            
            int myIndex = this->parent->getChildNumberContaining(this->value);
            this->parent->children[myIndex] = this; 
            this->fixHeightsRecursively();
        }

        #ifdef DEBUG_MODE
        template<typename T>
        void OctTree<T>::printHelper(const OctTreeNode *node, int indent) const
        {
            if(node == nullptr)
            {
                std::cout << "nullptr" << std::endl;
                return;
            }
            if(node->isLeaf)
            {
                std::cout << node->value << std::endl;
            }
            else
            {
                node->print();
            }
            int minNull = -1;
            for(int i = 0; i < CHILDREN - 1; i++)
            {
                if(node->children[i] != nullptr)
                {
                    if(minNull != -1)
                    {
                        for(int j = 0; j < indent; j++) std::cout << "\t";
                        if(minNull == i - 1)
                        {
                            std::cout << "[" << i-1 << "] nullptr" << std::endl;
                        }
                        else
                        {
                            std::cout << "[" << minNull << " - " << i-1 << "] nullptr" << std::endl;
                        }
                        minNull = -1;
                    }
                    for(int j = 0; j < indent; j++) std::cout << "\t";
                    std::cout << "[" << i << "] ";
                    this->printHelper(node->children[i], indent + 1);
                }
                else
                {
                    if(minNull == -1) minNull = i;
                }
            }
            if(minNull == -1)
            {
                for(int j = 0; j < indent; j++) std::cout << "\t";
                std::cout << "[" << (CHILDREN - 1) << "] ";
                this->printHelper(node->children[CHILDREN - 1], indent + 1);
            }
            else
            {
                if(node->children[CHILDREN - 1] != nullptr)
                {
                    for(int j = 0; j < indent; j++) std::cout << "\t";
                    std::cout << "[" << minNull << " - " << (CHILDREN - 2) << "] nullptr" << std::endl;
                    for(int j = 0; j < indent; j++) std::cout << "\t";
                    std::cout << "[" << (CHILDREN - 1) << "] ";
                    this->printHelper(node->children[CHILDREN - 1], indent + 1);
                }
                else
                {
                    for(int j = 0; j < indent; j++) std::cout << "\t";
                    std::cout << "[" << minNull << " - " << (CHILDREN - 1) << "] nullptr" << std::endl;
                }
            }
        }
        #endif // DEBUG_MODE

        template<typename T>
        template<typename U>
        typename OctTree<T>::OctTreeNode *OctTree<T>::tryInsert(const U &point)
        {
            if(this->getRoot() == nullptr)
            {
                throw MadVoro::Exception::MadVoroException("OctTree: root is null, the tree should be first configured with LL and UR");
            }
            if(!this->getRoot()->boundingBox.contains(point))
            {
                MadVoro::Exception::MadVoroException eo("OctTree: Value is outside the bounding box of the OctTree");
                eo.addEntry("Point", point);
                throw eo;
            }

            OctTreeNode *current = this->getRoot();
            while(current != nullptr)
            {
                // if we reached a leaf with the value `v`, start splitting until `v` and `point` are not in the same rectangle
                while(current->isLeaf)
                {
                    if(current->value == point)
                    {
                        return current;
                    }
                    current->splitNode();
                    current = current->parent;
                    current->id = (this->nodesNumber++);

                    int childIndex = current->getChildNumberContaining(point);
                    if(current->children[childIndex] == nullptr)
                    {
                        break;
                    }
                    current = current->children[childIndex];
                }
                // otherwise, determine the direction to go

                int childIndex = current->getChildNumberContaining(point);
                if(current->children[childIndex] == nullptr)
                {
                    return current->addLeafChild(childIndex, point);
                }
                current = current->children[childIndex];
            }
            MadVoro::Exception::MadVoroException eo("OctTree: A point could not be inserted to the OctTree");
            eo.addEntry("Point", point);
            throw eo;
        }

        template<typename T>
        typename OctTree<T>::OctTreeNode *OctTree<T>::removeLeaf(OctTreeNode *node)
        {
            // first, ensure it's indeed a leaf
            if(not node->isLeaf)
            {
                MadVoro::Exception::MadVoroException eo("OctTree remove: node could not be deleted, since its not a left");
                eo.addEntry("value", node->value);
                throw eo;
            }
            this->nodesNumber--;
            OctTreeNode *newParent = nullptr;
            if(node->parent == nullptr)
            {
                // remove the root
                this->root = nullptr;
                newParent = nullptr;
            }
            else
            {
                OctTreeNode *parent = node->parent;
                int myIndex;
                int numOfChildren;
                for(int i = 0; i < CHILDREN; i++)
                {
                    if(parent->children[i] != nullptr)
                    {
                        numOfChildren++;
                    }
                    if(parent->children[i] == node)
                    {
                        myIndex = i;
                        break;
                    }
                }
                if(numOfChildren == 0)
                {
                    throw MadVoro::Exception::MadVoroException("OctTree remove: should not reach here");
                }
                parent->children[myIndex] = nullptr;
                numOfChildren -= 1;
                if(numOfChildren == 0)
                {
                    // now a left
                    parent->isLeaf = true;
                }
                if(parent->isLeaf)
                {
                    newParent = this->removeLeaf(parent);
                }
                else
                {
                    newParent = parent;
                }
            }
            delete node;
            return newParent;
        }

        template<typename T>
        template<typename U, typename FilterFunction>
        std::vector<T> OctTree<T>::range(const Sphere<U> &sphere, size_t N, const FilterFunction &filter) const
        {
            std::vector<T> result;
            size_t resultSize = 0;
            this->nodes_stack.push_back(this->getRoot());

            while((not this->nodes_stack.empty()) and (resultSize < N))
            {
                const OctTreeNode *node = this->nodes_stack.back();
                this->nodes_stack.pop_back();

                if(node == nullptr)
                {
                    continue;
                }

                // DO NOT CHANGE THIS LINE TO "if `node->value` is in `sphere`"
                // that is because a leaf does not necessarily have to be a point (it can be a box, as in `DistributedOctTree`)
                if(not SphereBoxIntersection(node->boundingBox, sphere))
                {
                    continue;
                }
                
                if(node->isLeaf)
                {
                    if(filter(node->value))
                    {
                        result.push_back(node->value);
                        resultSize++;
                    }
                }
                else
                {
                    for(int i = 0; i < CHILDREN; i++)
                    {
                        this->nodes_stack.push_back(node->children[i]); // recursively iterate
                    }
                }
            }

            this->nodes_stack.clear();

            return result;
        }

        template<typename T>
        template<typename U, typename FilterFunction>
        std::pair<T, typename T::coord_type> OctTree<T>::getClosestPointInSphere(const Sphere<U> &sphere, const T &point, const FilterFunction &filter) const
        {
            this->nodes_stack.push_back(this->getRoot());

            T closestPoint;
            typename T::coord_type closestDistance = std::numeric_limits<typename T::coord_type>::max();

            while(not this->nodes_stack.empty())
            {
                const OctTreeNode *node = this->nodes_stack.back();
                this->nodes_stack.pop_back();

                if(node == nullptr)
                {
                    continue;
                }
                // calculate distance squared
                typename T::coord_type dist = node->boundingBox.distanceSquared(point);
                if((dist >= closestDistance) or (not SphereBoxIntersection(node->boundingBox, sphere)))
                {
                    continue;
                }
                // there might be a closer point in the subtrees
                if(node->isLeaf)
                {
                    if(filter(node->value))
                    {
                        // should not be ignored
                        closestPoint = node->value;
                        closestDistance = dist;
                    }
                }
                else
                {
                    for(int i = 0; i < CHILDREN; i++)
                    {
                        this->nodes_stack.push_back(node->children[i]);
                    }
                }
            }
            return {closestPoint, closestDistance};
        }

        template<typename T>
        const typename OctTree<T>::OctTreeNode *OctTree<T>::getNodeByDirections(const direction_t *directions) const
        {
            if(directions == nullptr)
            {
                return this->getRoot();
            }
            const OctTreeNode *current = this->getRoot();
            size_t i = 0;
            while(directions[i] != PATH_END_DIRECTION)
            {
                if(current == nullptr)
                {
                    break;
                };
                current = current->children[directions[i]];
                i++;
            }

            if(current == nullptr)
            {
                current = this->getRoot();
                MadVoro::Exception::MadVoroException eo("OctTree::getNodeByDirections: could not find the node");
                i = 0;
                while(directions[i] != PATH_END_DIRECTION)
                {

                    eo.addEntry("dir" + std::to_string(i), directions[i]);
                    if(current == nullptr)
                    {
                        eo.addEntry("Depth of null", i);
                        break;
                    }
                    else
                    {
                        eo.addEntry("node" + std::to_string(i), current->value);
                    }
                    current = current->children[directions[i]];
                    i++;
                }
                throw eo;
            }
            return current;
        }

        template<typename T>
        void OctTree<T>::getAllDecendantsHelper(const OctTreeNode *node, std::vector<T> &result) const
        {
            if(node == nullptr)
            {
                return;
            }
            if(node->isLeaf)
            {
                result.push_back(node->value);
            }
            for(int i = 0; i < CHILDREN; i++)
            {
                this->getAllDecendantsHelper(node->children[i], result);
            }
        }

        #define EPSILON 1e-12

        template<typename T>
        template<typename U>
        std::pair<T, typename T::coord_type> OctTree<T>::getClosestPointInfo(const U &point, bool includeSelf) const
        {
            this->nodes_stack.push_back(this->getRoot());

            T closestPoint;
            typename T::coord_type closestDistance = std::numeric_limits<typename T::coord_type>::max();
            
            while(!this->nodes_stack.empty())
            {
                const OctTreeNode *node = this->nodes_stack.back();
                this->nodes_stack.pop_back();

                if(node == nullptr)
                {
                    continue;
                }
                const T &closestPointInBox = node->boundingBox.closestPoint(point);
                // calculate distance squared
                typename T::coord_type dist = 0;
                for(int i = 0; i < DIM; i++)
                {
                    dist += (closestPointInBox[i] - point[i]) * (closestPointInBox[i] - point[i]);
                }
                if(dist >= closestDistance)
                {
                    continue;
                }
                // there might be a closer point in the subtrees
                if(node->isLeaf)
                {
                    if(not includeSelf and node->value == point)
                    {
                        // don't check that point (otherwise the distance is 0...)
                        continue;
                    }
                    closestPoint = node->value;
                    closestDistance = dist;
                }
                else
                {
                    for(int i = 0; i < CHILDREN; i++)
                    {
                        this->nodes_stack.push_back(node->children[i]);
                    }
                }
            }
            return {closestPoint, sqrt(closestDistance)};
        }
    }
}

template<typename T>
std::vector<T> MadVoro::DataStructure::OctTree<T>::getAllDecendants(const OctTreeNode *node) const
{
    std::vector<T> result;
    this->getAllDecendantsHelper(node, result);
    return result;
};

template<typename T>
template<typename InputIterator>
MadVoro::DataStructure::OctTree<T>::OctTree(const T &ll, const T &ur, const InputIterator &first, const InputIterator &last): OctTree(ll, ur)
{
    for(InputIterator it = first; it != last; it++)
    {
        this->insert(*it);
    }
};

template<typename T>
template<typename U>
T &MadVoro::DataStructure::OctTree<T>::findValue(const U &point)
{
    OctTreeNode *node = this->tryFind(point);
    if(node == nullptr)
    {
        MadVoro::Exception::MadVoroException eo("OctTree: could not find a point");
        eo.addEntry("Value", point);
        throw eo;
    }
    return node->value;
}

template<typename T>
template<typename U>
typename MadVoro::DataStructure::OctTree<T>::OctTreeNode *MadVoro::DataStructure::OctTree<T>::findNode(const U &point)
{
    OctTreeNode *node = this->tryFind(point);
    if(node == nullptr)
    {
        MadVoro::Exception::MadVoroException eo("OctTree: could not find a point");
        eo.addEntry("Value", point);
        throw eo;
    }
    return node;
}

#endif // _OCTTREE_HPP

