#ifndef HILBERT_RECTANGULAR_CONVERTOR_3D_HPP
#define HILBERT_RECTANGULAR_CONVERTOR_3D_HPP

#include "../HilbertConvertor3D.hpp"

namespace MadVoro
{
    namespace DataStructure
    {
        template<int max_leaf_ranks>        
        class HilbertTree3D;
    }
}

namespace MadVoro
{
    /**
     * see here the algorithm: https://github.com/jakubcerveny/gilbert
    */
    class HilbertRectangularConvertor3D : public HilbertConvertor3D
    {
        template<int max_leaf_ranks>        
        friend class MadVoro::DataStructure::HilbertTree3D;

    private:
        using direction_t = long int;

        struct DirectionPoint3D
        {
            direction_t x, y, z;

            #ifdef DEBUG_MODE
            friend std::ostream &operator<<(std::ostream &os, const DirectionPoint3D &args)
            {
                return os << "(" << args.x << ", " << args.y << ", " << args.z << ")";
            }
            #endif // DEBUG_MODE
        };

        struct RecursionArguments
        {
            DirectionPoint3D startPoint;
            DirectionPoint3D a;
            DirectionPoint3D b;
            DirectionPoint3D c;

            #ifdef DEBUG_MODE
            friend std::ostream &operator<<(std::ostream &os, const RecursionArguments &args)
            {
                return os << "startPoint = " << args.startPoint << ", a = " << args.a << ", b = " << args.b << ", c = " << args.c;
            }
            #endif // DEBUG_MODE
        };

        Point3D step;
        Geometry::BoundingBox<Point3D> spaceBoundingBox;
        DirectionPoint3D div;
        hilbert_index_t total_points_num;

    public:
        explicit HilbertRectangularConvertor3D(const Point3D &ll, const Point3D &ur, size_t order);
        
        inline hilbert_index_t getHilbertSize() const{return this->total_points_num;};
        
        void changeOrder(size_t order) override;
        
        hilbert_index_t xyz2d(coord_t x, coord_t y, coord_t z) const override;
            
        Point3D d2xyz(hilbert_index_t d) const override;
            
    private:
        std::vector<RecursionArguments> getRecursionArguments(const RecursionArguments &args) const;
        bool d2xyz_helper(const RecursionArguments &args, hilbert_index_t requested_d, hilbert_index_t &current_d, Point3D &result) const;
        bool xyz2d_helper_base(const DirectionPoint3D &startPoint, size_t steps, const DirectionPoint3D &direction, const DirectionPoint3D &requested_point, hilbert_index_t &current_d) const;
        bool xyz2d_helper(const RecursionArguments &args, const DirectionPoint3D &requested_point, hilbert_index_t &current_d) const;
        std::pair<DirectionPoint3D, DirectionPoint3D> getBoundingBox(const RecursionArguments &args) const;
        Point3D WidthHeightDepthToXYZ(direction_t width, direction_t height, direction_t depth) const;
    };
}

#endif // HILBERT_RECTANGULAR_CONVERTOR_3D_HPP