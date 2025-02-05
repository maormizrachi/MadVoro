#ifndef HILBERT_CONVERTOR_3D_HPP
#define HILBERT_CONVERTOR_3D_HPP

#ifdef DEBUG_MODE
    #include <iostream>
#endif // DEBUG_MODE
#include "elementary/Point3D.hpp" // for Point3D
#include "ds/utils/geometry.hpp" // for BoundingBox<Point3D>
#include "hilbertTypes.h"

#define MAX_HILBERT_ORDER 19

namespace MadVoro
{
    class HilbertConvertor3D
    {
    protected:
        using coord_t = Point3D::coord_type; // coordinate type

        Point3D ll, ur;
        size_t order;

    public:
        explicit HilbertConvertor3D(const Point3D &ll, const Point3D &ur, size_t order);
        
        virtual void changeOrder(size_t order) = 0;
        
        virtual hilbert_index_t xyz2d(coord_t x, coord_t y, coord_t z) const = 0;
        
        virtual inline hilbert_index_t xyz2d(const Point3D &point) const{return this->xyz2d(point.x, point.y, point.z);};
        
        virtual Point3D d2xyz(hilbert_index_t d) const = 0;
        
        inline size_t getOrder() const{return this->order;};
    };

    inline HilbertConvertor3D::HilbertConvertor3D(const Point3D &ll, const Point3D &ur, size_t order)
        : ll(ll), ur(ur), order(order)
    {}
}

#endif // HILBERT_CONVERTOR_3D_HPP