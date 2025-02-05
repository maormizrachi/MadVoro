#ifndef HILBERT_ORDINARY_CONVERTOR_3D_HPP
#define HILBERT_ORDINARY_CONVERTOR_3D_HPP

#include "../HilbertConvertor3D.hpp"
#include "../HilbertOrder3D.hpp"

namespace MadVoro
{
    class HilbertOrdinaryConvertor3D : public HilbertConvertor3D
    {
    public:
        explicit HilbertOrdinaryConvertor3D(const Point3D &ll, const Point3D &ur, size_t order);
        
        void changeOrder(size_t order) override
        {
            this->order = order;
        }
        
        hilbert_index_t xyz2d(coord_t x, coord_t y, coord_t z) const override;
            
        Point3D d2xyz(hilbert_index_t d) const override;
        
    private:
        HilbertCurve3D curve;
        Point3D length;
    };

    inline HilbertOrdinaryConvertor3D::HilbertOrdinaryConvertor3D(const Point3D &ll, const Point3D &ur, size_t order)
        : HilbertConvertor3D(ll, ur, order)
    {
        this->length = this->ur - this->ll; 
    }

    inline Point3D HilbertOrdinaryConvertor3D::d2xyz(hilbert_index_t d) const 
    {
        throw MadVoro::Exception::MadVoroException("HilbertOrdinaryConvertor3D::d2xyz: not implemented yet");
    }

    inline hilbert_index_t HilbertOrdinaryConvertor3D::xyz2d(coord_t x, coord_t y, coord_t z) const
    {
        Point3D vec;
        vec.x = (x - this->ll[0]) / this->length[0];
        vec.y = (y - this->ll[1]) / this->length[1];
        vec.z = (z - this->ll[2]) / this->length[2];
        return this->curve.Hilbert3D_xyz2d(vec, this->order);
    }
}

#endif // HILBERT_ORDINARY_CONVERTOR_3D_HPP