/*! \file HilbertOrder3D.hpp
\brief Hilbert 3D-space filling curve
\author Itai Linial
*/

// Need t do, add option to define number of points the same for early break, add protection against two points the same

#ifndef HILBERTORDER3D_HPP
#define HILBERTORDER3D_HPP 1
#include <vector>
#include <boost/array.hpp>
#include <algorithm>
#include <ctime>
#include <iostream>
#include "elementary/Point3D.hpp"
#include "HilbertOrder3D_Utils.hpp"

#define NUMBER_OF_SHAPES 24
#define MAX_ROTATION_LENGTH 5
#define PI 3.14159
#define MAX_HILBERT_RECURSIVE_CALLS 250

namespace MadVoro
{
	//! \brief The elementary Hilbert Curve shape
	class HilbertCurve3D_shape
	{
	public:
		//! \brief Constructor
		HilbertCurve3D_shape();
		//! \brief Comparison:
	//! \param shape Other shape
	//! \return Result of comparison
		bool operator==(const HilbertCurve3D_shape & shape) const;
		//! \brief An array of the 7 unit vector steps defining the shape:
		boost::array<Point3D, 7> m_vShapePoints;

	};

	//! \brief The Hilbert Curve class:
	class HilbertCurve3D
	{
	public:
		//! \brief Constructor
		HilbertCurve3D(void);
	//! \brief Calculate the Hilbert curve distance of a given point, given a required number of iterations:
	//! \param rvPoint Point
	//! \param numOfIterations Number of iterations
	//! \return Position on curve
		unsigned long long int Hilbert3D_xyz2d(Point3D const & rvPoint, int numOfIterations) const;

	private:
		// Rotate a shape according to a given rotation scheme (in-place):
		void RotateShape(int iShapeIndex, std::vector<int> vAxes);
		// Rotate a shape according to rotation index, and return the rotated shape:
		void RotateShape(HilbertCurve3D_shape const & roShape, HilbertCurve3D_shape & roShapeOut, int iRotationIndex);
		/*!
		\brief Returns the rotation scheme, according to a rotation index
		\param piRotation - a pointer to the output rotation scheme vector
		\param iRotationIndex - the desired rotation index
		\return The rotation scheme length (the size of the array given by piRotation)
		*/
		int GetRotation(int * piRotation, int iRotationIndex);
		// Find the index of a given shape object:
		int FindShapeIndex(const HilbertCurve3D_shape & roShape);
		// Create the recursion rule:
		void BuildRecursionRule();
		// Create the shape order, for all shapes (the order of octants):
		void BuildShapeOrder();

		// Stores all rotated shapes:
		boost::array<HilbertCurve3D_shape, NUMBER_OF_SHAPES> m_vRotatedShapes;
		// Stores all rotation schemes:
		boost::array < vector<int>, NUMBER_OF_SHAPES > m_vRotations;

		// An array of the 8 integers defining the recursion rule of the shape:
		boost::array< boost::array<int, 8> , NUMBER_OF_SHAPES> m_vShapeRecursion;

		// A 2x2x2 matrix indicating the 3 dimensional shape order
		// array< array<int , 8 > , NUMBER_OF_SHAPES > m_mShapeOrder;
		int m_mShapeOrder[NUMBER_OF_SHAPES][2][2][2];
	};

	/*!
	\brief Returns the 3D-Hilbert curve ordering
	\param cor The points
	\return The Hilbert order indices
	*/
	std::vector<std::size_t> HilbertOrder3D(std::vector<Point3D> const& cor);

	/*! \brief Get the reorderd list of indices
	\param cor List of points
	\param ll Low left corner
	\param ur Upper right corner
	\param Hmax Upper index of real indices
	\return List of indices
	*/
	std::vector<std::size_t> GetGlobalHibertIndeces(std::vector<Point3D> const& cor,Point3D const& ll, Point3D const& ur,size_t &Hmax);
}
#endif // HILBERTORDER_HPP
