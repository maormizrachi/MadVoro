#ifndef Point3D_HPP
#define Point3D_HPP 1

#include <vector>
#include <limits>
#include <cmath>
#include <cmath>
#include <math.h>
#include "utils/utils.hpp"

#ifdef MADVORO_WITH_MPI
	#include "mpi/serialize/Serializer.hpp"
#endif // MADVORO_WITH_MPI

using std::vector;

#define SIGN(x) ((x > 0) - (x < 0))

#ifndef EPSILON
	#define EPSILON 1e-12
#endif // EPSILON

namespace MadVoro
{
	//! \brief 3D Mathematical vector
	class Point3D
				#ifdef MADVORO_WITH_MPI
					: public MadVoro::MPI::Serializable
				#endif // MADVORO_WITH_MPI
	{
	public:
		using coord_type = double;

		//! \brief Component in the x direction
		double x;

		//! \brief Component in the y direction
		double y;

		//! \brief Component in the z direction
		double z;

		/*! \brief Class constructor
		\param ix x Component
		\param iy y Component
		\param iz z Component
		*/
		Point3D(double ix, double iy, double iz);

		/*! \brief Null constructor
		\details Sets all components to 0
		*/
		explicit Point3D(void);

		/*! \brief Class copy constructor
		\param other Other vector
		*/
		template<typename VectorType>
		Point3D(const VectorType &other): Point3D(other[0], other[1], other[2]){}

		~Point3D(void) = default;

		/*! \brief Set vector components
		\param ix x Component
		\param iy y Component
		\param iz z Component
		*/
		void Set(double ix, double iy, double iz);

	/*! \brief Indexed access to member
		\param index Member index
		\return Reference to member
	*/
		double& operator[](size_t index);

	/*! \brief Indexed access to member
		\param index Member index
		\return Value of member
	*/
		double operator[](size_t index) const;

		/*! \brief Addition
		\param v Vector to be added
		\return Reference to sum
		*/
		Point3D& operator+=(Point3D const& v);

		/*! \brief Subtraction
		\param v Vector to be subtracted
		\return Difference
		*/
		Point3D& operator-=(Point3D const& v);

		/*! \brief Assignment operator
		\param v Vector to be copied
		\return The assigned value
		*/
		template<typename VectorType>
		Point3D& operator=(const VectorType& v);

		/*! \brief Scalar product
		\param s Scalar
		\return Reference to the vector multiplied by scalar
		*/
		Point3D& operator*=(double s);

		/*! \brief Compare 3D-Vectors (up to an arbitrary precision)
		\param v Vector to be compared to
		\return True/False - according to the comparison results.
		*/
		bool operator==(Point3D const& v) const;

		bool operator!=(Point3D const& v) const;

		/*! \brief Rotates the vector around the X axes
		\param a Angle of rotation (in radians)
		*/
		void RotateX(double a);

		/*! \brief Rotates the vector around the Y axes
		\param a Angle of rotation (in radians)
		*/
		void RotateY(double a);

		/*! \brief Rotates the vector around the Z axes
		\param a Angle of rotation (in radians)
		*/
		void RotateZ(double a);

		/*! \brief Integer round of the vector's entries
		*/
		void Round();

		#ifdef MADVORO_WITH_MPI
			size_t dump(MadVoro::MPI::Serializer *serializer) const override;

			size_t load(const MadVoro::MPI::Serializer *serializer, size_t byteOffset) override;
		#endif // MADVORO_WITH_MPI

		friend std::ostream &operator<<(std::ostream &stream, const Point3D &vec);

		friend std::istream &operator>>(std::istream &stream, Point3D &vec);

		static const Point3D max(void);

		static const Point3D min(void);
	};

	/*! \brief Assignment operator
	\param v Vector to be copied
	\return The assigned value
	*/
	template<typename VectorType>
	#ifdef __INTEL_COMPILER
	#pragma omp declare simd
	#endif
	inline Point3D& Point3D::operator=(const VectorType& v)
	{
		x = v[0];
		y = v[1];
		z = v[2];
		return *this;
	}

	/*! \brief Norm of a vector
	\param v Three dimensional vector
	\return Norm of v
	*/
	#ifdef __INTEL_COMPILER
	#pragma omp declare simd
	#endif
	double abs(Point3D const& v);

	/*! \brief Norm of a vector, less accurate
	\param v Three dimensional vector
	\return Norm of v
	*/
	double fastabs(Point3D const& v);

	/*! \brief Term by term addition
	\param v1 First vector
	\param v2 Second vector
	\return Sum
	*/
	Point3D operator+(Point3D const& v1, Point3D const& v2);

	/*! \brief Term by term subtraction
	\param v1 First vector
	\param v2 Second vector
	\return Difference
	*/
	Point3D operator-(Point3D const& v1, Point3D const& v2);

	/*! \brief Scalar product
	\param d Scalar
	\param v Vector
	\return Three dimensional vector
	*/
	Point3D operator*(double d, Point3D const& v);

	/*! \brief Scalar product
	\param v Vector
	\param d Scalar
	\return Three dimensional vector
	*/
	Point3D operator*(Point3D const& v, double d);

	/*! \brief Scalar division
	\param v Vector
	\param d Scalar
	\return Three dimensional vector
	*/
	Point3D operator/(Point3D const& v, double d);

	/*! \brief Scalar product of two vectors
	\param v1 3D vector
	\param v2 3D vector
	\return Scalar product of v1 and v2
	*/
	double ScalarProd(Point3D const& v1, Point3D const& v2);

	/*! \brief Returns the angle between two vectors (in radians)
	\param v1 First vector
	\param v2 Second vector
	\return Angle (radians)
	*/
	double CalcAngle(Point3D const& v1, Point3D const& v2);

	/*! \brief Calculates the projection of one vector in the direction of the second
	\param v1 First vector
	\param v2 Direction of the projection
	\return Component of v1 in the direction of v2
	*/
	double Projection(Point3D const& v1, Point3D const& v2);

	/*! \brief Rotates a 3D-vector around the X axis
	\param v Vector
	\param a  (in radians)
	\return Rotated vector
	*/
	Point3D RotateX(Point3D const& v, double a);

	/*! \brief Rotates a 3D-vector around the Y axis
	\param v Vector
	\param a  (in radians)
	\return Rotated vector
	*/
	Point3D RotateY(Point3D const& v, double a);

	/*! \brief Rotates a 3D-vector around the Z axis
	\param v Vector
	\param a  (in radians)
	\return Rotated vector
	*/
	Point3D RotateZ(Point3D const& v, double a);

	/*! \brief Reflect vector
	\param v Vector
	\param normal Normal to the reflection plane
	\return Reflection of v about axis
	*/
	Point3D Reflect(Point3D const& v, Point3D const& normal);

	/*! \brief Calculates the distance between two vectors
	\param v1 First vector
	\param v2 Second vector
	\return distance between v1 and v2
	*/
	double distance(Point3D const& v1, Point3D const& v2);

	/*! \brief Returns the cross product of two vectors
	\param v1 First vector
	\param v2 Second vector
	\return Cross product between v1 and v2
	*/
	Point3D CrossProduct(Point3D const& v1, Point3D const& v2);

	/*! \brief Cross product
	\param v1 First vector
	\param v2 Second vector
	\param res result
	*/
	void CrossProduct(Point3D const& v1, Point3D const& v2,Point3D &res);

	/*! \brief Normalise vector
	\param vec Vector
	\return Normalised vector
	*/
	Point3D normalize(Point3D const& vec);

	/*! \brief Splits a vector of 3D points to components
	\param vIn Input vector of 3D points
	\param vX Vector of x coordinates (out)
	\param vY Vector of y coordinates (out)
	\param vZ Vector of z coordinates (out)
	*/
	void Split(vector<Point3D> const & vIn, vector<double> & vX, vector<double> & vY, vector<double> & vZ);

	template<>
	#ifdef __INTEL_COMPILER
	#pragma omp declare simd
	#endif
	inline Point3D::Point3D(const Point3D &v): Point3D(v.x, v.y, v.z)
	{}

	/*! \brief Assignment operator
	\param v Vector to be copied
	\return The assigned value
	*/
	template<>
	#ifdef __INTEL_COMPILER
	#pragma omp declare simd
	#endif
	inline Point3D& Point3D::operator=<Point3D>(const Point3D& v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}

	template<>
	inline Point3D::Point3D(const double &x): Point3D(x, x, x){}
}

#endif // Point3D_HPP

