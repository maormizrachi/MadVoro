#ifndef MADVORO_EXAMPLE_VECTOR3D_HPP
#define MADVORO_EXAMPLE_VECTOR3D_HPP

#include <iostream>

struct Vector3D
{
    using coord_type = double;

    double x, y, z;

    Vector3D() : x(0), y(0), z(0) {}
    Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}

    // Required by spatial_ds (BoundingBox, OctTree, Sphere)
    double &operator[](size_t i) { return (&x)[i]; }
    double operator[](size_t i) const { return (&x)[i]; }

    // Required by spatial_ds and MeshDecomposer3D
    friend bool operator==(const Vector3D &a, const Vector3D &b) { return a.x == b.x && a.y == b.y && a.z == b.z; }
    friend Vector3D operator+(const Vector3D &a, const Vector3D &b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
    friend Vector3D operator-(const Vector3D &a, const Vector3D &b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
    friend Vector3D operator*(const Vector3D &v, double s) { return {v.x * s, v.y * s, v.z * s}; }
    friend Vector3D operator*(double s, const Vector3D &v) { return {v.x * s, v.y * s, v.z * s}; }
    friend Vector3D operator/(const Vector3D &v, double s) { return {v.x / s, v.y / s, v.z / s}; }
    friend Vector3D &operator+=(Vector3D &a, const Vector3D &b) { a.x += b.x; a.y += b.y; a.z += b.z; return a; }
    friend Vector3D &operator-=(Vector3D &a, const Vector3D &b) { a.x -= b.x; a.y -= b.y; a.z -= b.z; return a; }
    friend Vector3D &operator*=(Vector3D &v, double s) { v.x *= s; v.y *= s; v.z *= s; return v; }
    friend Vector3D &operator/=(Vector3D &v, double s) { v.x /= s; v.y /= s; v.z /= s; return v; }

    friend std::ostream &operator<<(std::ostream &os, const Vector3D &v)
    {
        return os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    }
};

#endif // MADVORO_EXAMPLE_VECTOR3D_HPP
