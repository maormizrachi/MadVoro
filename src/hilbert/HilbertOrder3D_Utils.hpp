/*! \file HilbertOrder3D_Utils.hpp
\brief Hilbert 3D Order - Utility functions
\author Itai Linial
*/

#ifndef HILBERTORDER3D_UTILS_HPP
#define HILBERTORDER3D_UTILS_HPP 1

#include "elementary/Point3D.hpp"
#include <vector>
#include <algorithm>
#include <limits>
#include "hilbertTypes.h"
#include "utils/utils.hpp"

namespace MadVoro
{
	/*!
	\brief Estimate the number of iterations required in the Hilbert Curve, according to the number of points
	\param cor The points
	\return The estimated number of required iterations
	*/
	int EstimateHilbertIterationNum(std::vector<Point3D> const& cor);

	/*!
	\brief Scale a vector of 3D points to the unit-cube
	\param vPointsIn The input points
	\param vPointsOut (out) The output points
	*/
	void AdjustPoints(std::vector<Point3D> const & vPointsIn, std::vector<Point3D> & vPointsOut);

	/*!
	\brief Find points with same indeces
	\param vD_sorted The input points, sorted
	\param vOut (out) The list of points with same indeces
	*/
	void FindEqualIndices(std::vector<unsigned long long int> const & vD_sorted, std::vector<std::vector<std::size_t> > & vOut);

	/*! \brief Indices order after sorting of the values vector
	\param values Vector to be sorted
	\param indices Sorted list of indices
	*/
	template <typename T>
	void ordered(std::vector<T> const& values, std::vector<std::size_t> &indices)
	{
		indices.resize(values.size());
		MadVoro::Utils::sort_index(values,indices);
	}

	/*! \brief Reorder a vector according to an index vector (obtained from the 'ordered' function)
	\param v Vector to be sorted
	\param order Sorted list of indices
	*/
	template< class T >
	void reorder(std::vector<T> & v, std::vector<std::size_t> const & order)
	{
		std::vector<T> vCopy = v;
		for (std::size_t ii = 0; ii < v.size(); ++ii)
		{
			v[ii] = vCopy[order[ii]];
		}
		return;
	}
}

#endif // HILBERTORDER3D_UTILS_HPP
