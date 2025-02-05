#ifndef HILBERT_TREE_3D
#define HILBERT_TREE_3D

#define DEBUG_MODE

#ifdef DEBUG_MODE
    #include <iostream>
#endif // DEBUG_MODE

#include <vector>
#include <boost/container/flat_set.hpp>
#include <boost/container/small_vector.hpp>
#include <mpi.h>
#include "ds/utils/geometry.hpp"
#include "HilbertRectangularConvertor3D.hpp"

namespace MadVoro
{
    namespace DataStructure
    {
        template<typename T>
        using BoundingBox = MadVoro::Geometry::BoundingBox<T>;
        
        template<typename T>
        using Sphere = MadVoro::Geometry::Sphere<T>;

        #define DEFAULT_RANKS_IN_LEAVES 4
        #define UNDEFINED_OWNER -1

        template<int max_leaf_ranks = DEFAULT_RANKS_IN_LEAVES>
        class HilbertTree3DNode
        {
        public:
            explicit HilbertTree3DNode(HilbertTree3DNode *parent): parent(parent){};
            
            BoundingBox<Point3D> boundingBox;
            hilbert_index_t d_start, d_end;
            size_t num_points; // number of points in the subtree
            std::vector<HilbertTree3DNode*> children;
            bool is_leaf;
            boost::container::small_vector<int, max_leaf_ranks> owners;
            HilbertTree3DNode *parent;
        };

        template<int max_leaf_ranks = DEFAULT_RANKS_IN_LEAVES>
        class HilbertTree3D
        {
        public:
            using Node = HilbertTree3DNode<max_leaf_ranks>;

        private:
            using RanksSet = boost::container::flat_set<int>;

        private:
            Node *root;
            mutable std::vector<const Node*> nodes_stack;
            MPI_Comm comm;
            int rank, size;
            #ifdef DEBUG_MODE
                const HilbertConvertor3D *convertor;
            #endif // DEBUG_MODE

            void buildTreeHelper(Node *currentNode, const typename HilbertRectangularConvertor3D::RecursionArguments &current_args, hilbert_index_t &current_d, const HilbertRectangularConvertor3D *convertor, const std::vector<hilbert_index_t> &responsibilityRange);
            
            void buildTree(const HilbertRectangularConvertor3D *convertor, const std::vector<hilbert_index_t> &responsibilityRange);

            #ifdef DEBUG_MODE
                void printHelper(const Node *node, int tabs = 0) const;
            #endif // DEBUG_MODE

        public:
            HilbertTree3D(const HilbertRectangularConvertor3D *convertor, const std::vector<hilbert_index_t> &responsibilityRange, const MPI_Comm &comm = MPI_COMM_WORLD): comm(comm)
            {
                MPI_Comm_rank(this->comm, &this->rank);
                MPI_Comm_size(this->comm, &this->size);
                #ifdef DEBUG_MODE
                    this->convertor = convertor;
                #endif // DEBUG_MODE
                this->buildTree(convertor, responsibilityRange);
            }

            ~HilbertTree3D();

            template<typename U>
            RanksSet getIntersectingRanks(const Sphere<U> &sphere) const;

            template<typename U>
            inline RanksSet getIntersectingRanks(const U &point, typename U::coord_type radius) const
            {
                return this->getIntersectingRanks(Sphere<U>(point, radius));
            }

            template<typename U>
            std::vector<std::pair<typename Point3D::coord_type, typename Point3D::coord_type>> getClosestFurthestPointsByRanks(const U &point) const;

            #ifdef DEBUG_MODE
                void print() const{this->printHelper(this->root);};
            #endif // DEBUG_MODE

            std::vector<BoundingBox<Point3D>> getRankBoundingBoxes(int _rank) const;
            
            inline std::vector<BoundingBox<Point3D>> getMyBoundingBoxes() const
            {
                return this->getRankBoundingBoxes(this->rank);
            }

            std::vector<std::vector<BoundingBox<Point3D>>> getBoundingBoxesOfRanks(void) const;

            std::vector<const Node*> getValuesIf(const std::function<bool(const Node*)> ifOpenFunction, const std::function<bool(const Node*)> &ifAddValueFunction) const;
        };

        template<int max_ranks_per_leaf>
        HilbertTree3D<max_ranks_per_leaf>::~HilbertTree3D()
        {
            this->nodes_stack.push_back(this->root);
            while(not this->nodes_stack.empty())
            {
                const Node *node = this->nodes_stack.back();
                this->nodes_stack.pop_back();

                if(node == nullptr)
                {
                    continue;
                }

                if(!node->is_leaf)
                {
                    for(const Node *child : node->children)
                    {
                        this->nodes_stack.push_back(child);
                    }
                }

                delete node;
            }

            this->root = nullptr;
        }

        #ifdef DEBUG_MODE
        template<int max_ranks_per_leaf>
        void HilbertTree3D<max_ranks_per_leaf>::printHelper(const Node *node, int tabs) const
        {
            if(node == nullptr)
            {
                return;
            }
            for(int i = 0; i < tabs; i++) std::cout << "\t";
            std::cout << "LL = " << node->boundingBox.getLL() << ", UR = " << node->boundingBox.getUR() << ", d: " << node->d_start << " - " << node->d_end << " (points: " << node->num_points << "). ";

            size_t numOwners = node->owners.size();

            if(numOwners == 0)
            {
                std::cout << "No explicit owner";
            }
            else
            {
                if(numOwners == 1)
                {
                    std::cout << "Owner: " << node->owners[0];
                }
                else
                {
                    std::cout << "Owners: ";
                    for(size_t i = 0; i < numOwners - 1; i++)
                    {
                        std::cout << node->owners[i] << ", ";
                    }
                    std::cout << node->owners[numOwners - 1];
                }
            }

            std::cout << ((node->is_leaf)? " LEAF" : " NOT LEAF");

            std::cout << std::endl;

            for(const Node *child : node->children)
            {
                this->printHelper(child, tabs + 1);
            }
        }
        #endif // DEUBG_MODE

        template<int max_ranks_per_leaf>
        void HilbertTree3D<max_ranks_per_leaf>::buildTreeHelper(Node *currentNode, const typename HilbertRectangularConvertor3D::RecursionArguments &current_args, hilbert_index_t &current_d, const HilbertRectangularConvertor3D *convertor, const std::vector<hilbert_index_t> &responsibilityRange)
        {
            using DirectionPoint3D = HilbertRectangularConvertor3D::DirectionPoint3D;
            using RecursionArguments = HilbertRectangularConvertor3D::RecursionArguments;
            using direction_t = HilbertRectangularConvertor3D::direction_t;

            if(currentNode == nullptr)
            {
                return;
            }

            const DirectionPoint3D &startPoint = current_args.startPoint;
            const DirectionPoint3D &a = current_args.a;
            const DirectionPoint3D &b = current_args.b;
            const DirectionPoint3D &c = current_args.c;
            direction_t width = std::abs(a.x + a.y + a.z);
            direction_t height = std::abs(b.x + b.y + b.z);
            direction_t depth = std::abs(c.x + c.y + c.z);

            size_t num_points = width * height * depth;

            currentNode->num_points = num_points;
            
            currentNode->d_start = current_d;
            currentNode->d_end = current_d + num_points;
            
            // set bounding box
            const std::pair<DirectionPoint3D, DirectionPoint3D> &boundingBox = convertor->getBoundingBox(current_args);
            const DirectionPoint3D &ll = boundingBox.first;       
            const DirectionPoint3D &ur = boundingBox.second;       
            currentNode->boundingBox = BoundingBox<Point3D>(convertor->WidthHeightDepthToXYZ(ll.x, ll.y, ll.z), convertor->WidthHeightDepthToXYZ(ur.x, ur.y, ur.z));

            std::pair<int, int> ranksMatching = {0, this->size - 1};


            if(current_d >= responsibilityRange.back())
            {
                ranksMatching = {this->size-1, this->size-1};
            }
            else
            {
                for(int index = 0; index < this->size; index++)
                {
                    if(currentNode->d_start <= responsibilityRange[index])
                    {
                        ranksMatching.first = index;
                        break;
                    }
                }

                for(int index = ranksMatching.first; index < this->size; index++)
                {
                    if((currentNode->d_end - 1) <= responsibilityRange[index])
                    {
                        ranksMatching.second = index;
                        break;
                    }

                }
            }

            Point3D newLL = std::numeric_limits<typename Point3D::coord_type>::max() * Point3D(1, 1, 1);
            Point3D newUR = std::numeric_limits<typename Point3D::coord_type>::lowest() * Point3D(1, 1, 1);

            if((ranksMatching.second - ranksMatching.first) < max_ranks_per_leaf)
            {
                // don't have to call recursively
                currentNode->is_leaf = true;
                size_t numOwners = ranksMatching.second - ranksMatching.first + 1;
                currentNode->owners.resize(numOwners);
                for(size_t i = 0; i < numOwners; i++)
                {
                    currentNode->owners[i] = ranksMatching.first + i;
                }
                current_d += num_points; // one big step for d
            }
            else
            {
                currentNode->is_leaf = false;
                currentNode->owners.resize(0);

                // should call recursively
                direction_t dax = SIGN(a.x), day = SIGN(a.y), daz = SIGN(a.z);
                direction_t dbx = SIGN(b.x), dby = SIGN(b.y), dbz = SIGN(b.z);
                direction_t dcx = SIGN(c.x), dcy = SIGN(c.y), dcz = SIGN(c.z);

                // for base case:
                bool baseCase = false;
                DirectionPoint3D baseCaseUnitDirection;
                direction_t baseCaseLength;
                // check for base cases
                if(height == 1 and depth == 1)
                {
                    baseCase = true;
                    baseCaseUnitDirection = {dax, day, daz};
                    baseCaseLength = width;
                }

                if(width == 1 and depth == 1)
                {
                    baseCase = true;
                    baseCaseUnitDirection = {dbx, dby, dbz};
                    baseCaseLength = height;
                }

                if(width == 1 and height == 1)
                {
                    baseCase = true;
                    baseCaseUnitDirection = {dcx, dcy, dcz};
                    baseCaseLength = depth;
                }

                if(baseCase)
                {
                    direction_t x = startPoint.x, y = startPoint.y, z = startPoint.z;
                    for(int i = 0; i < baseCaseLength; i++)
                    {
                        currentNode->children.push_back(new Node(currentNode));
                        this->buildTreeHelper(currentNode->children.back(), {{x, y, z}, {dax, day, daz}, {dbx, dby, dbz}, {dcx, dcy, dcz}}, current_d, convertor, responsibilityRange);
                        x += baseCaseUnitDirection.x;
                        y += baseCaseUnitDirection.y;
                        z += baseCaseUnitDirection.z;
                    }
                }
                else
                {
                    for(const RecursionArguments &nextArgs : convertor->getRecursionArguments(current_args))
                    {
                        currentNode->children.push_back(new Node(currentNode));
                        this->buildTreeHelper(currentNode->children.back(), nextArgs, current_d, convertor, responsibilityRange);
                    }      

                    for(const Node *child : currentNode->children)
                    {
                        const Point3D &childLL = child->boundingBox.getLL();
                        const Point3D &childUR = child->boundingBox.getUR();
                        
                        for(int j = 0; j < DIM; j++)
                        {
                            newLL[j] = std::min<typename Point3D::coord_type>(newLL[j], childLL[j]);
                            newUR[j] = std::max<typename Point3D::coord_type>(newUR[j], childUR[j]);
                        }
                    }    
                    newLL -= Point3D(EPSILON, EPSILON, EPSILON);
                    newUR += Point3D(EPSILON, EPSILON, EPSILON);
                    currentNode->boundingBox.setBounds(newLL, newUR);
                }
            }
        }

        template<int max_ranks_per_leaf>
        void HilbertTree3D<max_ranks_per_leaf>::buildTree(const HilbertRectangularConvertor3D *convertor, const std::vector<hilbert_index_t> &responsibilityRange)
        {
            using DirectionPoint3D = HilbertRectangularConvertor3D::DirectionPoint3D;
            using RecursionArguments = HilbertRectangularConvertor3D::RecursionArguments;
            using direction_t = HilbertRectangularConvertor3D::direction_t;
            hilbert_index_t d = 0;

            this->root = new Node(nullptr);
            this->buildTreeHelper(this->root, {{0, 0, 0}, {convertor->div.x, 0, 0}, {0, convertor->div.y, 0}, {0, 0, convertor->div.z}}, d, convertor, responsibilityRange);

            if(d != convertor->total_points_num)
            {
                MadVoro::Exception::MadVoroException eo("HilbertTree3D::buildTree: d != convertor->total_points_num");
                eo.addEntry("d", d);
                eo.addEntry("total points num", convertor->total_points_num);
                throw eo;
            }
        }

        template<int max_ranks_per_leaf>
        template<typename U>
        typename HilbertTree3D<max_ranks_per_leaf>::RanksSet HilbertTree3D<max_ranks_per_leaf>::getIntersectingRanks(const Sphere<U> &sphere) const
        {
            RanksSet result;
            this->nodes_stack.push_back(this->root);

            while(not this->nodes_stack.empty())
            {
                const Node *node = this->nodes_stack.back();
                this->nodes_stack.pop_back();

                if(node == nullptr)
                {
                    continue;
                }

                if(not SphereBoxIntersection(node->boundingBox, sphere))
                {
                    continue;
                }

                if(node->is_leaf)
                {
                    size_t numOwners = node->owners.size();
                    for(size_t i = 0; i < numOwners; i++)
                    {
                        result.insert(node->owners[i]);
                    }
                }
                else
                {
                    for(const Node *child : node->children)
                    {
                        this->nodes_stack.push_back(child); // recursively iterate
                    }
                }
            }

            if(result.empty())
            {
                MadVoro::Exception::MadVoroException eo("In HilbertTree3D::getIntersectingRanks: result is empty, (should at least contain the rank itself)");
                eo.addEntry("Rank", rank);
                eo.addEntry("Sphere", sphere);
                eo.addEntry("Root Bounding Box", this->root->boundingBox);
                #ifdef DEBUG_MODE
                    eo.addEntry("d of sphere center", this->convertor->xyz2d(sphere.center));
                    this->print();
                #endif // DEBUG_MODE
                throw eo;
            }

            return result;
        }

        template<int max_ranks_per_leaf>
        template<typename U>
        std::vector<std::pair<typename Point3D::coord_type, typename Point3D::coord_type>> HilbertTree3D<max_ranks_per_leaf>::getClosestFurthestPointsByRanks(const U &point) const
        {
            using coord_type = typename Point3D::coord_type;

            const coord_type &maxVal = std::numeric_limits<coord_type>::max();
            const coord_type &minVal = std::numeric_limits<coord_type>::lowest();
            
            std::pair<Point3D, Point3D> initialPair = std::make_pair<Point3D, Point3D>(Point3D(maxVal, maxVal, maxVal), Point3D(minVal, minVal, minVal));
            std::vector<std::pair<coord_type, coord_type>> distances(this->size, {maxVal, minVal});

            this->nodes_stack.push_back(this->root);

            Point3D closestPoint, furthestPoint;
            while(not this->nodes_stack.empty())
            {
                const Node *node = this->nodes_stack.back();
                this->nodes_stack.pop_back();

                if(node == nullptr)
                {
                    continue;
                }
                if(!node->is_leaf)
                {
                    for(const Node *child : node->children)
                    {
                        this->nodes_stack.push_back(child);
                    }
                    continue;
                }
                // node is a value node
                closestPoint = node->boundingBox.closestPoint(point);
                furthestPoint = node->boundingBox.furthestPoint(point);
                coord_type closestDist = 0, furthestDist = 0;
                for(int i = 0; i < DIM; i++)
                {
                    closestDist += (closestPoint[i] - point[i]) * (closestPoint[i] - point[i]);
                    furthestDist += (furthestPoint[i] - point[i]) * (furthestPoint[i] - point[i]);
                }

                size_t numOwners = node->owners.size();
                for(size_t i = 0; i < numOwners; i++)
                {
                    int owner = node->owners[i];
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

        template<int max_leaf_ranks>
        std::vector<BoundingBox<Point3D>> HilbertTree3D<max_leaf_ranks>::getRankBoundingBoxes(int _rank) const
        {
            std::vector<BoundingBox<Point3D>> result;

            this->nodes_stack.push_back(this->root);
            while(not this->nodes_stack.empty())
            {
                const Node *node = this->nodes_stack.back();
                this->nodes_stack.pop_back();
                if(node == nullptr)
                {
                    continue;
                }
                if(node->is_leaf)
                {
                    if(std::find(node->owners.begin(), node->owners.end(), _rank) != node->owners.end())
                    {
                        result.push_back(node->boundingBox);
                    }
                }
                else
                {
                    for(const Node *child : node->children)
                    {
                        this->nodes_stack.push_back(child);
                    }
                }
            }
            return result;
        }

        template<int max_leaf_ranks>
        std::vector<std::vector<BoundingBox<Point3D>>> HilbertTree3D<max_leaf_ranks>::getBoundingBoxesOfRanks(void) const
        {
            std::vector<std::vector<BoundingBox<Point3D>>> result(this->size);

            this->nodes_stack.push_back(this->root);
            while(not this->nodes_stack.empty())
            {
                const Node *node = this->nodes_stack.back();
                this->nodes_stack.pop_back();
                if(node == nullptr)
                {
                    continue;
                }
                if(node->is_leaf)
                {
                    for(const int &_rank : node->owners)
                    {
                        result[_rank].push_back(node->boundingBox);
                    }
                }
                else
                {
                    for(const Node *child : node->children)
                    {
                        this->nodes_stack.push_back(child);
                    }
                }
            }
            return result;
        }

        template<int max_leaf_ranks>
        std::vector<const typename HilbertTree3D<max_leaf_ranks>::Node*> HilbertTree3D<max_leaf_ranks>::getValuesIf(const std::function<bool(const Node*)> ifOpenFunction, const std::function<bool(const Node*)> &ifAddValueFunction) const
        {
            std::vector<const Node*> nodes = {this->root};
            nodes.reserve(this->getDepth() * max_leaf_ranks);

            std::vector<Node*> result;

            while(not nodes.empty())
            {
                const Node *node = nodes.back();
                nodes->pop_back();
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
                    for(const Node *child : node->children)
                    {
                        nodes.push_back(child);
                    }
                }
            }
            return result;
        }
    }
}

#endif // HILBERT_TREE_3D