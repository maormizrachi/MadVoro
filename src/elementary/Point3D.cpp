#include "Point3D.hpp"

namespace
{
	static double my_round(double val)
	{    
		return floor(val + 0.5);
	}
}

MadVoro::Point3D::Point3D(double ix, double iy, double iz): x(ix), y(iy), z(iz) {};

MadVoro::Point3D::Point3D(void): Point3D(0, 0, 0){};

void MadVoro::Point3D::Set(double ix, double iy, double iz) 
{
    x = ix;
    y = iy;
    z = iz;
}

double& MadVoro::Point3D::operator[](size_t index)
{
    double *values[3] = {&this->x, &this->y, &this->z};
    return *values[index];
}

/*! \brief Indexed access to member
\param index Member index
\return Value of member
*/
double MadVoro::Point3D::operator[](size_t index) const
{
    double values[3] = {this->x, this->y, this->z};
    return values[index];
}

#ifdef __INTEL_COMPILER
#pragma omp declare simd
#endif
MadVoro::Point3D& MadVoro::Point3D::operator+=(Point3D const& v)
{
    x += v.x;
    y += v.y;
    z += v.z;
    return *this;
}

#ifdef __INTEL_COMPILER
#pragma omp declare simd
#endif
MadVoro::Point3D& MadVoro::Point3D::operator-=(Point3D const& v)
{
    x -= v.x;
    y -= v.y;
    z -= v.z;
    return *this;
}

MadVoro::Point3D& MadVoro::Point3D::operator*=(double s)
{
    x *= s;
    y *= s;
    z *= s;
    return *this;
}

bool MadVoro::Point3D::operator==(Point3D const& v) const
{
    // Note - since working with double precision, two vectors are assumed to be "equal",
    // if their coordinates agree up to precision EPSILON
    return (std::abs(x - v.x) < EPSILON) && (std::abs(y - v.y) < EPSILON) && (std::abs(z - v.z) < EPSILON);
}

bool MadVoro::Point3D::operator!=(Point3D const& v) const
{
    return !this->operator==(v);
}

void MadVoro::Point3D::RotateX(double a)
{
    Point3D v;
    v.x = x;
    v.y = y*cos(a) - z*sin(a);
    v.z = y*sin(a) + z*cos(a);

    *this = v;
}

void MadVoro::Point3D::RotateY(double a)
{
    Point3D v;
    v.x = x*cos(a) + z*sin(a);
    v.y = y;
    v.z = -x*sin(a) + z*cos(a);

    *this = v;
}

void MadVoro::Point3D::RotateZ(double a)
{
    Point3D v;
    v.x = x*cos(a) - y*sin(a);
    v.y = x*sin(a) + y*cos(a);
    v.z = z;

    *this = v;
}

void MadVoro::Point3D::Round()
{
    x = my_round(x);
    y = my_round(y);
    z = my_round(z);
}

#ifdef MADVORO_WITH_MPI
    size_t MadVoro::Point3D::dump(MadVoro::MPI::Serializer *serializer) const
    {
        size_t bytes = 0;
        bytes += serializer->insert(this->x);
        bytes += serializer->insert(this->y);
        bytes += serializer->insert(this->z);
        return bytes;
    }

    size_t MadVoro::Point3D::load(const MadVoro::MPI::Serializer *serializer, size_t byteOffset)
    {
        size_t bytes = 0;
        bytes += serializer->extract(this->x, byteOffset);
        bytes += serializer->extract(this->y, byteOffset + bytes);
        bytes += serializer->extract(this->z, byteOffset + bytes);
        return bytes;
    }
#endif // MADVORO_WITH_MPI

static const MadVoro::Point3D max(void)
{
    return MadVoro::Point3D(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
}

static const MadVoro::Point3D min(void)
{
    return MadVoro::Point3D(std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest(), std::numeric_limits<double>::lowest());
}

double MadVoro::abs(Point3D const& v)
{
    return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

double MadVoro::fastabs(Point3D const& v)
{
    return MadVoro::Utils::fastsqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

MadVoro::Point3D MadVoro::operator+(Point3D const& v1, Point3D const& v2)
{
    Point3D res;
    res.x = v1.x + v2.x;
    res.y = v1.y + v2.y;
    res.z = v1.z + v2.z;
    return res;
}

MadVoro::Point3D MadVoro::operator-(Point3D const& v1, Point3D const& v2)
{
    Point3D res;
    res.x = v1.x - v2.x;
    res.y = v1.y - v2.y;
    res.z = v1.z - v2.z;
    return res;
}

MadVoro::Point3D MadVoro::operator*(double d, Point3D const& v)
{
    Point3D res;
    res.x = v.x * d;
    res.y = v.y * d;
    res.z = v.z * d;
    return res;
}

MadVoro::Point3D MadVoro::operator*(Point3D const& v, double d)
{
    return d*v;
}

MadVoro::Point3D MadVoro::operator/(Point3D const& v, double d)
{
    Point3D res;
    res.x = v.x / d;
    res.y = v.y / d;
    res.z = v.z / d;
    return res;
}

double MadVoro::ScalarProd(Point3D const& v1, Point3D const& v2)
{
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

double MadVoro::CalcAngle(Point3D const& v1, Point3D const& v2)
{
    return std::acos(std::max(-1.0, std::min(1.0, ScalarProd(v1, v2) / (abs(v1) * abs(v2)))));
}

double MadVoro::Projection(Point3D const& v1, Point3D const& v2)
{
    return ScalarProd(v1, v2) / abs(v2);
}

MadVoro::Point3D MadVoro::RotateX(Point3D const& v, double a)
{
    Point3D res;
    res.x = v.x;
    res.y = v.y*cos(a) - v.z*sin(a);
    res.z = v.y*sin(a) + v.z*cos(a);
    return res;
}

MadVoro::Point3D MadVoro::RotateY(Point3D const& v, double a)
{
    Point3D res;
    res.x = v.x*cos(a) + v.z*sin(a);
    res.y = v.y;
    res.z = -v.x*sin(a) + v.z*cos(a);
    return res;
}

MadVoro::Point3D MadVoro::RotateZ(Point3D const& v, double a)
{
    Point3D res;
    res.x = v.x*cos(a) - v.y*sin(a);
    res.y = v.x*sin(a) + v.y*cos(a);
    res.z = v.z;
    return res;
}

MadVoro::Point3D MadVoro::Reflect(Point3D const& v, Point3D const& normal)
{
    return v - 2 * ScalarProd(v, normal)*normal / ScalarProd(normal,normal);
}

double MadVoro::distance(Point3D const& v1, Point3D const& v2)
{
    return abs(v1 - v2);
}

MadVoro::Point3D MadVoro::CrossProduct(Point3D const& v1, Point3D const& v2)
{
    return Point3D(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}

void MadVoro::CrossProduct(Point3D const& v1, Point3D const& v2,Point3D &res)
{
    res.x = v1.y*v2.z - v1.z*v2.y;
    res.y = v1.z*v2.x - v1.x*v2.z;
    res.z = v1.x*v2.y - v1.y*v2.x;
}

MadVoro::Point3D MadVoro::normalize(Point3D const& vec)
{
    double l = abs(vec);
    return vec / l;
}

void MadVoro::Split(const std::vector<Point3D> &vIn, std::vector<double> &vX, std::vector<double> &vY, std::vector<double> &vZ)
{
    vX.resize(vIn.size());
    vY.resize(vIn.size());
    vZ.resize(vIn.size());

    for (std::size_t ii = 0; ii < vIn.size(); ++ii)
    {
        vX[ii] = vIn[ii].x;
        vY[ii] = vIn[ii].y;
        vZ[ii] = vIn[ii].z;
    }
    return;
}

namespace MadVoro
{
    std::ostream &operator<<(std::ostream &stream, const Point3D &vec)
    {
        stream << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")";
        return stream;
    }

    std::istream &operator>>(std::istream &stream, Point3D &vec)
    {
        std::string str;
        std::getline(stream, str, '(');
        std::getline(stream, str, ',');
        vec.x = std::stod(str);
        std::getline(stream, str, ',');
        vec.y = std::stod(str);
        std::getline(stream, str, ')');
        vec.z = std::stod(str);
        return stream;
    }
}