#include "Vector3D.hpp"

namespace MadVoro
{
    Vector3D operator+(Vector3D const& v1, Vector3D const& v2)
    {
        Vector3D res;
        res.x = v1.x + v2.x;
        res.y = v1.y + v2.y;
        res.z = v1.z + v2.z;
        return res;
    }

    Vector3D operator-(Vector3D const& v1, Vector3D const& v2)
    {
        Vector3D res;
        res.x = v1.x - v2.x;
        res.y = v1.y - v2.y;
        res.z = v1.z - v2.z;
        return res;
    }

    Vector3D operator*(double d, Vector3D const& v)
    {
        Vector3D res;
        res.x = v.x * d;
        res.y = v.y * d;
        res.z = v.z * d;
        return res;
    }

    Vector3D operator*(Vector3D const& v, double d)
    {
        return d*v;
    }

    Vector3D operator/(Vector3D const& v, double d)
    {
        Vector3D res;
        res.x = v.x / d;
        res.y = v.y / d;
        res.z = v.z / d;
        return res;
    }

    std::ostream &operator<<(std::ostream& os, const Vector3D& v)
    {
        return os << "(" << v.x << ", " << v.y << ", " << v.z <<  ")";
    }
}
