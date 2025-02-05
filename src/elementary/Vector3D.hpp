#ifndef MADVORO_POINT_HPP
#define MADVORO_POINT_HPP

#include <iostream>

namespace MadVoro
{
    class Vector3D
    {
    public:
        double x, y, z;

        Vector3D(double x, double y, double z): x(x), y(y), z(z) {}
        
        Vector3D(): Vector3D(0, 0, 0) {}

        ~Vector3D() = default;

    };

    /*! \brief Term by term addition
    \param v1 First vector
    \param v2 Second vector
    \return Sum
    */
    Vector3D operator+(Vector3D const& v1, Vector3D const& v2);

    /*! \brief Term by term subtraction
    \param v1 First vector
    \param v2 Second vector
    \return Difference
    */
    Vector3D operator-(Vector3D const& v1, Vector3D const& v2);

    /*! \brief Scalar product
    \param d Scalar
    \param v Vector
    \return Three dimensional vector
    */
    Vector3D operator*(double d, Vector3D const& v);

    /*! \brief Scalar product
    \param v Vector
    \param d Scalar
    \return Three dimensional vector
    */
    Vector3D operator*(Vector3D const& v, double d);

    /*! \brief Scalar division
    \param v Vector
    \param d Scalar
    \return Three dimensional vector
    */
    Vector3D operator/(Vector3D const& v, double d);
    
    std::ostream &operator<<(std::ostream& os, const Vector3D& v);
};

#endif // MADVORO_POINT_HPP