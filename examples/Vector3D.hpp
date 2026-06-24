#ifndef EXAMPLE_VECTOR3D_HPP
#define EXAMPLE_VECTOR3D_HPP

#include <iostream>
#include <cmath>

/**
 * @brief A simple 3D vector type for use with MadVoro examples.
 *
 * Your own project should provide its own Vec3 type with public
 * `double x, y, z` members and that is constructible via {x, y, z}.
 */
struct Vector3D
{
    double x, y, z;

    Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}
    Vector3D() : x(0), y(0), z(0) {}

    double& operator[](size_t i) { return (&x)[i]; }
    double  operator[](size_t i) const { return (&x)[i]; }

    Vector3D& operator+=(const Vector3D &o) { x+=o.x; y+=o.y; z+=o.z; return *this; }
    Vector3D& operator-=(const Vector3D &o) { x-=o.x; y-=o.y; z-=o.z; return *this; }
    Vector3D& operator*=(double d) { x*=d; y*=d; z*=d; return *this; }
    Vector3D& operator/=(double d) { x/=d; y/=d; z/=d; return *this; }

    bool operator==(const Vector3D &o) const { return x==o.x && y==o.y && z==o.z; }
    bool operator!=(const Vector3D &o) const { return !(*this == o); }

    friend Vector3D operator+(const Vector3D &a, const Vector3D &b) { return {a.x+b.x, a.y+b.y, a.z+b.z}; }
    friend Vector3D operator-(const Vector3D &a, const Vector3D &b) { return {a.x-b.x, a.y-b.y, a.z-b.z}; }
    friend Vector3D operator*(double d, const Vector3D &v) { return {d*v.x, d*v.y, d*v.z}; }
    friend Vector3D operator*(const Vector3D &v, double d) { return d*v; }
    friend Vector3D operator/(const Vector3D &v, double d) { return {v.x/d, v.y/d, v.z/d}; }

    friend std::ostream &operator<<(std::ostream &os, const Vector3D &v)
    {
        return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    }
};

inline double ScalarProduct(const Vector3D &a, const Vector3D &b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

inline Vector3D CrossProduct(const Vector3D &a, const Vector3D &b)
{
    return {a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z, a.x*b.y - a.y*b.x};
}

inline double abs(const Vector3D &v)
{
    return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

inline Vector3D Normalize(const Vector3D &v)
{
    double len = abs(v);
    return (len > 0) ? v / len : v;
}

#endif // EXAMPLE_VECTOR3D_HPP
