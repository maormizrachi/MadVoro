#ifndef _KDTREE_HPP
#define _KDTREE_HPP

#include <vector>
#include <functional>
#include <utility>

#include "ds/utils/geometry.hpp"

#ifdef KD_DEBUG_MODE
#include <iostream>
#endif // KD_DEBUG_MODE

namespace MadVoro
{
    namespace DataStructure
    {
        template<typename T>
        using BoundingBox = MadVoro::Geometry::BoundingBox<T>;
        
        template<typename T>
        using Sphere = MadVoro::Geometry::Sphere<T>;

        /**
         * A kd tree class.
         * @param T the type of point
         * @param D the dimension
        */
        template<typename T, int D>
        class KDTree
        {
            using DefaultFilterFunction = std::function<bool(const T&)>;

        private:
            class KDTreeNode
            {
            friend class KDTree;

            public:

                KDTreeNode(const T &value, KDTreeNode *parent); // for non-root nodes
                KDTreeNode(const T &value, const T &ll, const T &ur); // for root

                T value;
                KDTreeNode *left, *right;
                KDTreeNode *parent;
                bool isValue; // whether a real value in the tree, or merely an auxiliary node
                int partitionAxis;
                int height;
                int depth;
                BoundingBox<T> boundingBox;

            private:
                enum _ChildSide
                {
                    LEFT,
                    RIGHT
                };

                void updateHeightRecursively();
                inline _ChildSide getChildSide(const T &value) const
                {
                    if(value[this->partitionAxis] <= this->value[this->partitionAxis])
                    {
                        return _ChildSide::LEFT;
                    }
                    return _ChildSide::RIGHT;
                } 
            };

        public:
            explicit inline KDTree(const T &ll, const T &ur): root(nullptr), size(0), boundingBox(BoundingBox<T>(ll, ur)){};
            inline ~KDTree(){this->deleteTree();};
            
            template<typename U, typename FilterFunction = DefaultFilterFunction>
            std::vector<T> range(const Sphere<U> &sphere, size_t N = std::numeric_limits<size_t>::max(), const FilterFunction &filter = [](const T&){return true;}) const;

            template<typename U, typename FilterFunction = DefaultFilterFunction>
            std::pair<T, typename T::coord_type> getClosestPointInSphere(const Sphere<U> &sphere, const T &point, const FilterFunction &filter = [](const T&){return true;}) const;

            inline bool find(const T &value) const{return this->tryFind(value) != nullptr;};
            inline bool insert(const T &value){return this->tryInsert(value) != nullptr;};
            inline void deleteTree(){this->deleteHelper(this->getRoot()); this->setRoot(nullptr);};
            inline size_t getSize() const{return this->size;};
            #ifdef KD_DEBUG_MODE
            void print() const{this->printHelper(this->getRoot(), 0);};
            #endif // KD_DEBUG_MODE

        private:
            KDTreeNode *root;
            size_t size;
            BoundingBox<T> boundingBox;

            const KDTreeNode *tryFind(const T &point) const;
            inline KDTreeNode *tryFind(const T &point){return const_cast<KDTreeNode*>(std::as_const(*this).tryFind(point));};
            const KDTreeNode *tryFindParent(const T &point) const;
            inline KDTreeNode *tryFindParent(const T &point){return const_cast<KDTreeNode*>(std::as_const(*this).tryFindParent(point));};
            KDTreeNode *tryInsert(const T &value);
            void deleteHelper(KDTreeNode *root);

            inline KDTreeNode *getRoot(){return this->root;};
            inline const KDTreeNode *getRoot() const{return this->root;};
            inline void setRoot(KDTreeNode *newRoot){this->root = newRoot;};
            #ifdef KD_DEBUG_MODE
            void printHelper(const KDTreeNode *node, int tabs) const;
            #endif // KD_DEBUG_MODE
        };

        template<typename T, int D>
        void KDTree<T, D>::KDTreeNode::updateHeightRecursively()
        {
            KDTreeNode *current = this;
            while(current != nullptr)
            {
                int leftHeight = (current->left == nullptr)? -1 : current->left->height;
                int rightHeight = (current->right == nullptr)? -1 : current->right->height;
                current->height = std::max<int>(leftHeight, rightHeight) + 1;
                current = current->parent;
            }
        }

        template<typename T, int D>
        KDTree<T, D>::KDTreeNode::KDTreeNode(const T &value, const T &ll, const T &ur): value(value), left(nullptr), right(nullptr), parent(nullptr), isValue(true)
        {
            this->boundingBox = BoundingBox<T>(ll, ur);
            this->depth = 0;
            this->partitionAxis = 0;
        }

        template<typename T, int D>
        KDTree<T, D>::KDTreeNode::KDTreeNode(const T &value, KDTreeNode *parent): value(value), left(nullptr), right(nullptr), parent(parent), isValue(true)
        {
            if(parent == nullptr)
            {
                this->partitionAxis = 0;
                this->depth = 0;
            }
            else
            {
                this->partitionAxis = (this->parent->partitionAxis + 1) % DIM;
                this->depth = this->parent->depth + 1;

                // hang on parent correct size
                T my_ll = this->parent->boundingBox.getLL(), my_ur = this->parent->boundingBox.getUR();
            
                if(this->parent->getChildSide(value) == _ChildSide::LEFT)
                {
                    this->parent->left = this;
                    my_ur[this->parent->partitionAxis] = this->parent->value[this->parent->partitionAxis];
                }
                else
                {
                    this->parent->right = this;
                    my_ll[this->parent->partitionAxis] = this->parent->value[this->parent->partitionAxis];
                }
                this->boundingBox = BoundingBox<T>(my_ll, my_ur);
            }

            this->updateHeightRecursively();
        }

        template<typename T, int D>
        typename KDTree<T, D>::KDTreeNode *KDTree<T, D>::tryInsert(const T &value)
        {
            if(!this->boundingBox.contains(value))
            {
                // point is out of the tree range
                return nullptr;
            }
            KDTreeNode *parent = this->tryFindParent(value);
            KDTreeNode *node;
            if(parent == nullptr)
            {
                // root insertion
                node = new KDTreeNode(value, this->boundingBox.getLL(), this->boundingBox.getUR());
                this->setRoot(node);
            }
            else
            {
                node = new KDTreeNode(value, parent);
            }
            this->size++;
            return node;
        }

        template<typename T, int D>
        const typename KDTree<T, D>::KDTreeNode *KDTree<T, D>::tryFindParent(const T &value) const
        {
            const KDTreeNode *current = this->getRoot();
            while(current != nullptr)
            {
                if(current->getChildSide(value) == KDTreeNode::_ChildSide::LEFT)
                {
                    if(current->left == nullptr)
                    {
                        return current;
                    }
                    current = current->left; // go left
                }
                else
                {
                    if(current->right == nullptr)
                    {
                        return current;
                    }
                    current = current->right; // go right
                }
            }
            return nullptr;
        }

        template<typename T, int D>
        const typename KDTree<T, D>::KDTreeNode *KDTree<T, D>::tryFind(const T &value) const
        {
            const KDTreeNode *current = this->getRoot();
            while(current != nullptr)
            {
                if((current->value == value) and (current->isValue))
                {
                    return current;
                }

                if(current->getChildSide(value) == KDTreeNode::_ChildSide::LEFT)
                {
                    current = current->left; // go left
                }
                else
                {
                    current = current->right; // go right
                }
            }
            return nullptr;
        }

        template<typename T, int D>
        void KDTree<T, D>::deleteHelper(typename KDTree<T, D>::KDTreeNode *root)
        {
            if(root == nullptr)
            {
                return;
            }
            this->deleteHelper(root->left);
            this->deleteHelper(root->right);
            this->size--;
            delete root;
        }

        template<typename T, int D>
        template<typename U, typename FilterFunction>
        std::vector<T> KDTree<T, D>::range(const Sphere<U> &sphere, size_t N, const FilterFunction &filter) const
        {
            std::vector<T> result;
            size_t resultSize = 0;
            std::vector<const KDTreeNode*> nodes_stack;
            nodes_stack.push_back(this->getRoot());

            while((resultSize < N) and (not nodes_stack.empty()))
            {
                const KDTreeNode *node = nodes_stack.back();
                nodes_stack.pop_back();

                if(node == nullptr)
                {
                    continue;
                }

                if(not SphereBoxIntersection(node->boundingBox, sphere))
                {
                    continue;
                }

                if(node->isValue)
                {
                    if(filter(node->value) and sphere.contains(node->value))
                    {
                        result.push_back(node->value);
                        resultSize++;
                    }
                }
                nodes_stack.push_back(node->left);
                nodes_stack.push_back(node->right);
            }
            return result;
        }

        template<typename T, int D>
        template<typename U, typename FilterFunction>
        std::pair<T, typename T::coord_type> KDTree<T, D>::getClosestPointInSphere(const Sphere<U> &sphere, const T &point, const FilterFunction &filter) const
        {
            std::vector<const KDTreeNode*> nodes_stack;
            nodes_stack.push_back(this->getRoot());
            T closestPoint;
            typename T::coord_type closestDistance = std::numeric_limits<typename T::coord_type>::max();

            while(not nodes_stack.empty())
            {
                const KDTreeNode *node = nodes_stack.back();
                nodes_stack.pop_back();

                if(node == nullptr)
                {
                    continue;
                }
                T closestPointInBox = node->boundingBox.closestPoint(point);
                // calculate distance squared
                typename T::coord_type dist = 0;
                for(int i = 0; i < D; i++)
                {
                    dist += (closestPointInBox[i] - point[i]) * (closestPointInBox[i] - point[i]);
                }
                if((dist >= closestDistance) or (not SphereBoxIntersection(node->boundingBox, sphere)))
                {
                    continue;
                }
                // there might be a closer point in the subtrees
                if(node->isValue)
                {
                    if(filter(node->value) and sphere.contains(node->value))
                    {
                        // should not be ignored
                        closestPoint = node->value;
                        closestDistance = dist;
                    }
                }
                nodes_stack.push_back(node->left);
                nodes_stack.push_back(node->right);
            }
            return {closestPoint, closestDistance};
        }

        #ifdef KD_DEBUG_MODE
        template<typename T, int D>
        void KDTree<T, D>::printHelper(const KDTreeNode *node, int tabs) const
        {
            if(node == nullptr)
            {
                for(int i = 0; i < tabs; i++){std::cout << "\t";}
                std::cout << "nullptr" << std::endl;
                return;
            }
            for(int i = 0; i < tabs; i++){std::cout << "\t";}
            if(node->isValue)
            {
                std::cout << "\033[1;34m";
            }
            std::cout << node->value << " (height: " << node->height << ", height: " << node->height << ", [bounding box " << node->boundingBox.getLL() << " x " << node->boundingBox.getUR() << "])";
            std::cout << "\033[1;0m" << std::endl;
            for(int i = 0; i < tabs; i++){std::cout << "\t";}
            std::cout << "[L]" << std::endl;
            printHelper(node->left, tabs + 1);
            for(int i = 0; i < tabs; i++){std::cout << "\t";}
            std::cout << "[R]" << std::endl;
            printHelper(node->right, tabs + 1);
        }
        #endif // KD_DEBUG_MODE
    }
}

#endif // _KDTREE_HPP
