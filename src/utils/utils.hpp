/*!  \file utils.hpp
\brief Various useful functions
\author Almog Yalinewich
*/

#ifndef UTILS_HPP
#define UTILS_HPP 1

#include <vector>
#include <algorithm>
#include <cassert>
#include <iostream>
#include <boost/container/flat_map.hpp>
#include <limits>
#include <cmath>
#include "exception/SizeException.hpp"

namespace MadVoro
{
	namespace Utils
	{
		using std::vector;

		/*! \brief Checks if a number is numerically close to zero
		\param x Number
		\return True is the number is close enough to zero
		*/
		bool close2zero(const double x);

		/*! \brief Term by term addition for vectors
		\param v1 Right argument
		\param v2 Left argument
		\return The addition of the two vectors
		*/
		template <class T> vector<T> operator+
		(vector<T> const& v1, vector<T> const& v2);

		template <class T> vector<T> operator+
		(vector<T> const& v1, vector<T> const& v2)
		{
			assert(v1.size() == v2.size() && "Vectors must have the same length");
			vector<T> res(v1.size());
			for (std::size_t i = 0; i < v1.size(); ++i)
				res[i] = v1[i] + v2[i];
			return res;
		}

		/*! \brief Adds the same thing to all terms in a vector
		\param v Vector
		\param t Addition
		\return New vector with sums as terms
		*/
		template<class T> vector<T> operator+
		(const vector<T>& v, const T& t);

		template<class T> vector<T> operator+
		(const vector<T>& v, const T& t)
		{
			vector<T> res(v.size());
			for (std::size_t i = 0, endp = v.size(); i < endp; ++i)
				res[i] = v[i] + t;
			return res;
		}

		/*! \brief Adds the same thing to all terms in a vector
		\param v Vector
		\param t Addition
		\return New vector with sums as terms
		*/
		template<class T> vector<T> operator+
		(const T& t, const vector<T>& v);

		template<class T> vector<T> operator+
		(const T& t, const vector<T>& v)
		{
			return v + t;
		}

		/*! \brief Adds the same thing to all terms in a vector
		\param v Vector
		\param t Addition
		\return Reference to united vector
		*/
		template<class T> vector<T>& operator+=
		(vector<T>& v, const T& t);

		template<class T> vector<T>& operator+=
		(vector<T>& v, const T& t)
		{
			for (std::size_t i = 0, endp = v.size(); i < endp; ++i)
				v[i] += t;
			return v;
		}

		/*! \brief Term by term vector subtraction
		\param v1 Left argument
		\param v2 Right argument
		\return  The subtraction between the vectors
		*/
		template <class T> vector<T> operator-
		(vector<T> const& v1, vector<T> const& v2);

		template <class T> vector<T> operator-
		(vector<T> const& v1, vector<T> const& v2)
		{
			if (v1.size() != v2.size())
				throw MadVoro::Exception::MadVoroException("Vectors must have the same length");

			vector<T> res(v1.size());
			for (std::size_t i = 0; i < v1.size(); ++i)
				res[i] = v1[i] - v2[i];
			return res;
		}

		/*! \brief Multiplies all terms of a vector by a scalar
		\param d Scalar
		\param v vector
		\return The mulitplication of the vector
		*/
		template <class T> vector<T> operator*
		(double d, vector<T> const& v);

		template <class T> vector<T> operator*
		(double d, vector<T> const& v)
		{
			vector<T> res(v.size());
			for (std::size_t i = 0; i < v.size(); ++i)
				res[i] = d*v[i];
			return res;
		}
		/*!
		\brief BiLinear Interpolation
		\param x The x vector, assumed sorted
		\param y The y vector, assumed sorted
		\param xi The x interpolation location
		\param yi The y interpolation location
		\param data data(x,y) The function to interpolate
		\return f(xi, yi)
		*/
		template <typename T>
		T BiLinearInterpolation(std::vector<T> const& x, std::vector<T> const& y, std::vector<std::vector<T>> const& data, T const xi, T const yi);

		template <typename T>
		T BiLinearInterpolation(std::vector<T> const& x, std::vector<T> const& y, std::vector<std::vector<T>> const& data, T const xi, T const yi)
		{
			auto itx = std::upper_bound(x.begin(), x.end(), xi);
			if(itx == x.end())
			{
				if(xi > x.back())
					throw MadVoro::Exception::MadVoroException("X too large in BiLinearInterpolation");
				else
					--itx;
			}
			if(itx == x.begin())
				throw MadVoro::Exception::MadVoroException("X too small in BiLinearInterpolation");
			auto ity = std::upper_bound(y.begin(), y.end(), yi);
			if(ity == y.end())
			{
				if(yi > y.back())
					throw MadVoro::Exception::MadVoroException("Y too large in BiLinearInterpolation");
				else
					--ity;
			}
			if(ity == y.begin())
				throw MadVoro::Exception::MadVoroException("Y too small in BiLinearInterpolation");
			T const delta = 1.0 / ((*itx - *(itx - 1)) * (*ity - *(ity - 1)));
			size_t const idx = static_cast<size_t>(itx - x.begin());
			size_t const idy = static_cast<size_t>(ity - y.begin());
			return (data[idx - 1][idy - 1] * (*itx - xi) * (*ity - yi) + data[idx][idy - 1] * (xi - *(itx - 1)) * (*ity - yi)
				+ data[idx - 1][idy] * (*itx - xi) * (yi - *(ity -1)) + data[idx - 1][idy - 1] * (xi - *(itx - 1)) * (yi - *(ity -1))) * delta;
		}

		/*!
		\brief Linear Interpolation
		\param x The x vector, assumed sorted
		\param y y=f(x) vector
		\param xi The interpolation location
		\return f(xi)
		*/
		template <typename T>
		T LinearInterpolation(const vector<T> &x, const vector<T> &y, T xi);

		template <typename T>
		T LinearInterpolation(const vector<T> &x, const vector<T> &y, T xi)
		{
			typename vector<T>::const_iterator it = upper_bound(x.begin(), x.end(), xi);
			if (it == x.end())
			{
			//		if (x.back() == xi)
			if(MadVoro::Utils::close2zero(x.back()-xi))
					return y.back();
				std::cout << "X too large in LinearInterpolation, x_i " << xi << " max X " << x.back() << std::endl;
				throw;
			}
			if (it == x.begin())
			{
				if (*it < x.at(0))
				{
					std::cout << "X too small in LinearInterpolation, x_i " << xi << " min X " << x.at(0) << std::endl;
					throw;
				}
			}
			if (MadVoro::Utils::close2zero(*it - xi))
				return y[static_cast<std::size_t>(it - x.begin())];

			return y[static_cast<std::size_t>(it - x.begin())] + (xi - *it)*	(y[static_cast<std::size_t>(it - 1 - x.begin())] 
				- y[static_cast<std::size_t>(it - x.begin())]) / (*(it - 1) - *it);
		}

		/*! \brief Returns the minimal term in a vector
		\param v Vector
		\return The minimum of the vector
		*/
		double min(vector<double> const& v);

		/*! \brief returns the maximal term in a vector
		\param v Vector
		\return The maximum of the vector
		*/
		double max(vector<double> const& v);

		/*!
		\brief Removes the elements in v given by indeces
		\param v The vector to change
		\param indeces The cells to remove
		*/
		template <class T> void RemoveVector
		(vector<T> &v, vector<int> &indeces);

		template <class T> void RemoveVector
		(vector<T> &v, vector<int> &indeces)
		{
			if (indeces.empty())
				return;
			sort(indeces.begin(), indeces.end());
			vector<T> result;
			result.reserve(v.size() - indeces.size());
			int counter = 0;
			for (std::size_t i = 0; i < static_cast<std::size_t>(indeces.back()); ++i)
			{
				if (std::size_t(indeces[std::size_t(counter)]) == i)
					++counter;
				else
					result.push_back(v[i]);
			}
			for (std::size_t i = static_cast<std::size_t>(indeces.back()) + 1; i < v.size(); ++i)
				result.push_back(v[i]);
			v = result;
		}

		/*!
		\brief Removes the elements in v given by indeces
		\param v The vector to change
		\param indeces The cells to remove
		*/
		template <class T> void RemoveVector
		(vector<T> &v, vector<std::size_t> &indeces);

		template <class T> void RemoveVector
		(vector<T> &v, vector<std::size_t> &indeces)
		{
			if (indeces.empty())
				return;
			sort(indeces.begin(), indeces.end());
			vector<T> result;
			result.reserve(v.size() - indeces.size());
			std::size_t counter = 0;
			for (std::size_t i = 0; i < indeces.back(); ++i)
			{
				if (indeces[counter] == i)
					++counter;
				else
					result.push_back(v[i]);
			}
			for (std::size_t i = indeces.back() + 1; i < v.size(); ++i)
				result.push_back(v[i]);
			v = result;
		}
		/*!
		\brief Returns only the values with indeces in index
		\param v The vector to check
		\param index The indeces to return
		\return The reduced vector
		*/
		template<class T>
		vector<T> VectorValues(vector<T> const&v, vector<int> const &index)
		{
			if (index.empty() || v.empty())
				return vector<T>();

			vector<T> result(index.size());
			for (std::size_t i = 0; i < index.size(); ++i)
				result.at(i) = v.at(static_cast<std::size_t>(index.at(i)));
			return result;
		}

		/*! \brief Returns only the values with indeces in index
		\param v The vector to check
		\param index The indeces to return
		\return The reduced vector
		*/
		template <class T> vector<T> VectorValues
		(vector<T> const& v,
		vector<std::size_t> const &index)
		{
			if (index.empty() || v.empty())
				return vector<T>();

			vector<T> result(index.size());
			for (std::size_t i = 0; i < index.size(); ++i)
				result.at(i) = v.at(index.at(i));
			return result;
		}

		/*! \brief Reverses the vector
		\param v Vector
		*/
		template <class T> void FlipVector(vector<T> &v)
		{
			vector<T> temp(v);
			size_t N = temp.size();
			for (size_t i = 0; i < N; ++i)
				temp[i] = v[N - i - 1];
			v = temp;
		}

		/*!
		\brief Returns the sum of the vector
		\param v The vector to sum
		\return The sum
		*/
		template <class T> T VectorSum(vector<T> const&v);

		template <class T> T VectorSum(vector<T> const&v)
		{
			if (v.empty())
				return 0;
			T result = v[0];
			for (std::size_t i = 1; i < v.size(); ++i)
				result += v[i];
			return result;
		}

		/*!
		\brief Returns a vector containing only unique elements
		\param v The input vector, must be SORTED!!
		\return The unique vector
		*/
		template <class T> vector<T> unique(vector<T> const& v);

		template <class T> vector<T> unique(vector<T> const& v)
		{
			std::size_t n = v.size();
			vector<T> res;
			res.reserve(n);
			if (n == 0)
				return res;
			res.push_back(v[0]);
			for (typename vector<T>::const_iterator it = v.begin() + 1; it != v.end(); ++it)
				if (*it == *(it - 1))
					continue;
				else
					res.push_back(*it);
			return res;
		}

		/*!
		\brief Returns a vector containing only indeces of unique elements
		\param v The input vector, must be SORTED!!
		\return The unique vector indeces
		*/
		template <class T> vector<int> unique_index(vector<T> const& v);

		template <class T> vector<int> unique_index(vector<T> const& v)
		{
			if (v.empty())
				return vector<int>();

			vector<int> res;
			res.push_back(0);
			for (std::size_t i = 1; i < v.size(); ++i)
				if (v[i] != v[i - 1])
					res.push_back(static_cast<int>(i));
			return res;
		}

		/*!
		\brief Returns only elements from vector v which are not in vector list, assumes list is sorted
		\param v The vector to change
		\param list The values to be taken out of v
		\return The new vector
		*/
		template <class T> vector<T> RemoveList(vector<T> const&v, vector<T> const&list);

		template <class T> vector<T> RemoveList(vector<T> const&v, vector<T> const&list)
		{
			vector<T> res;
			for (std::size_t i = 0; i < v.size(); ++i)
				if (!std::binary_search(list.begin(), list.end(), v[i]))
					res.push_back(v[i]);
			return res;
		}

		/*!
		\brief Removes the first occurence of val inside a vector
		\param vec The vector to change
		\param val The value to remove
		*/
		template <class T> void RemoveVal(vector<T> &vec, T val);

		template <class T> void RemoveVal(vector<T> &vec, T val)
		{
			for (std::size_t i = 0; i < vec.size(); ++i)
			{
				if (vec[i] == val)
				{
					vec.erase(vec.begin() + static_cast<long>(i));
					return;
				}
			}
		}

		/*!
		\brief checks if val is in the vector
		\param vec The vector to check
		\param val The value to check inside the vector
		\return True if the value is in the vector, false otherwise
		*/
		template <class T> bool InVector(vector<T> const&vec, T val);

		template <class T> bool InVector(vector<T> const&vec, T val)
		{
			int n = vec.size();
			for (int i = 0; i < n; ++i)
				if (vec[i] == val)
					return true;
			return false;
		}

		/*!
		\brief Returns the index of val in the vector
		\param vec The vector to check
		\param val The value to look for
		\return The index of the vector which first equals to val, throws exception if not found
		*/
		template <class T> int IndexInVector(vector<T> const&vec, T val);

		template <class T> int IndexInVector(vector<T> const&vec, T val)
		{
			int n = vec.size();
			for (int i = 0; i < n; ++i)
				if (vec[i] == val)
					return i;
			throw MadVoro::Exception::MadVoroException("Value not found in vector");
		}
		/*!
		\brief Rearranges the vector according to the indeces
		\param v The vector to rearrange
		\param indeces The rearrangement indeces
		*/
		template <class T> void ReArrangeVector(vector<T> &v, vector<int> const& indeces);

		template <class T> void ReArrangeVector(vector<T> &v, vector<int> const& indeces)
		{
			const vector<T> temp = v;
			for (std::size_t i = 0; i < v.size(); ++i)
				v[i] = temp[static_cast<std::size_t>(indeces[i])];
		}

		/*!
		\brief Rearranges the vector according to the indeces
		\param v The vector to rearrange
		\param indeces The rearrangement indeces
		*/
		template <class T> void ReArrangeVector(vector<T> &v, vector<std::size_t> const& indeces);

		template <class T> void ReArrangeVector(vector<T> &v, vector<std::size_t> const& indeces)
		{
			const vector<T> temp = v;
			for (std::size_t i = 0; i < v.size(); ++i)
				v[i] = temp[indeces[i]];
		}

		namespace
		{
			template<class T> struct index_cmp
			{
				explicit index_cmp(const T _arr) : arr(_arr) {}
				bool operator()(const std::size_t a, const std::size_t b) const
				{
					return arr[a] < arr[b];
				}
			private:
				const T arr;
			};

			template<class T> struct Wrapper
			{
				index_cmp<vector<T> > *cmp;
				explicit Wrapper(index_cmp<vector<T> > *cmp2) : cmp(cmp2) {}
				bool operator() (const size_t lhs, const size_t rhs) const
				{
					return (*cmp)(lhs, rhs);
				}
			};
		}
		/*! \brief Returns the indeces of a sort
		\param arr The array to sort
		\param res The indeces of the sort that is given as the output
		*/
		template<class T> void sort_index(const vector<T> & arr, vector<int>& res);

		template<class T> void sort_index(const vector<T> & arr, vector<int>& res)
		{
			res.resize(arr.size());
			for (std::size_t i = 0; i < res.size(); ++i)
				res[i] = static_cast<int>(i);
			index_cmp<vector<T> > cmp(arr);
			Wrapper<T> w(&cmp);
			sort(res.begin(), res.end(), w);
		}

		/*! \brief Returns the indeces of a sort
		\param arr The array to sort
		\param res The indeces of the sort that is given as the output
		*/
		template<class T> void sort_index(const vector<T> & arr, vector<std::size_t>& res);

		template<class T> void sort_index(const vector<T> & arr, vector<std::size_t>& res)
		{
			res.resize(arr.size());
			for (std::size_t i = 0; i < res.size(); ++i)
				res[i] = i;
			index_cmp<vector<T> > cmp(arr);
			Wrapper<T> w(&cmp);
			sort(res.begin(), res.end(), w);
		}


		/*! \brief Returns the indeces of a sort
		\param arr The array to sort
		\return The indeces of the sort that is given as the output
		*/
		template<class T> vector<std::size_t> sort_index(const vector<T> & arr);

		template<class T> vector<std::size_t> sort_index(const vector<T> & arr)
		{
			if(arr.size() == 0)
				return vector<std::size_t>();
			vector<std::size_t> res(arr.size());
			for (std::size_t i = 0; i < res.size(); ++i)
				res[i] = i;
			index_cmp<vector<T> > cmp(arr);
			Wrapper<T> w(&cmp);
			sort(res.begin(), res.end(), w);
			return res;
		}

		namespace
		{
			template <class RAIter, class Compare>
			class PairComp
			{
			public:
				Compare comp;
			explicit PairComp(Compare comp_) : comp(comp_) {}
				bool operator() (const std::pair<std::size_t, RAIter>& a,
					const std::pair<std::size_t, RAIter>& b) const {
					return comp(*a.second, *b.second);
				}
			};
		}

		/*! \brief Returns the indeces of a sort
		\param iterBegin Starting iterator
		\param iterEnd End iterator
		\param comp The compare function
		\param indexes Output
		*/
		template <class RAIter, class Compare>
		void sort_index(RAIter iterBegin, RAIter iterEnd, Compare comp,
			std::vector<std::size_t>& indexes);

		template <class RAIter, class Compare>
		void sort_index(RAIter iterBegin, RAIter iterEnd, Compare comp,
			std::vector<std::size_t>& indexes)
		{

			std::vector< std::pair<std::size_t, RAIter> > pv;
			pv.reserve(iterEnd - iterBegin);

			RAIter iter;
			std::size_t k;
			for (iter = iterBegin, k = 0; iter != iterEnd; iter++, k++) {
				pv.push_back(std::pair<size_t, RAIter>(k, iter));
			}
			PairComp<RAIter, Compare> compy(comp);
			std::sort(pv.begin(), pv.end(), compy);

			indexes.resize(pv.size());
			for (std::size_t i = 0; i < pv.size(); ++i)
				indexes[i] = pv[i].first;
		}

		/*! \brief Concatenates two vectors
		\param v1 vector
		\param v2 vector
		\return vector
		*/
		template<class T> vector<T> join(vector<T> const& v1,
			vector<T> const& v2);

		template<class T> vector<T> join(vector<T> const& v1,
			vector<T> const& v2)
		{
			if (v1.empty() && v2.empty())
				return vector<T>();
			vector<T> res(v1.size() + v2.size());
			if (!v1.empty())
				copy(v1.begin(), v1.end(), res.begin());
			if (!v2.empty())
				copy(v2.begin(), v2.end(), res.begin() + static_cast<int>(v1.size()));
			return res;
		}

		//! \brief BinaryOperation Binary operator
		template<class T> class BinaryOperation
		{
		public:

			/*! \brief Evaluates the binary operation
			\param t1 First argument
			\param t2 Second argument
			\return Result
			*/
			virtual T operator()(T const& t1, T const& t2) const = 0;

			virtual ~BinaryOperation(void) {}
		};

		/*! \brief Applies a binary operation on every pair of values from two vectors
		\param v1 First vector
		\param v2 Second vector
		\param bin_op Binary operation
		\return Vector
		*/
		template<class T> vector<T> binary_unite(vector<T> const& v1,
			vector<T> const& v2,
			BinaryOperation<T> const& bin_op);

		template<class T> vector<T> binary_unite(vector<T> const& v1,
			vector<T> const& v2,
			BinaryOperation<T> const& bin_op)
		{
			assert(v1.size() == v2.size());

			vector<T> res(v1.size());
			for (std::size_t i = 0; i < v1.size(); ++i)
				res[i] = bin_op(v1[i], v2[i]);
			return res;
		}

		//! \brief Abstract class for an unary operation
		template<class T> class UnaryOperation;

		template<class T> class UnaryOperation
		{
		public:

			/*! \brief Evaluate operator
			\param t Argument
			\return Result of operation
			*/
			virtual T operator()(T const& t) const = 0;

			virtual ~UnaryOperation(void) {}
		};

		/*! \brief Applies an unary operator to all terms in a std::vector
		\param v Vector
		\param un_op Unary operator
		\return Vector
		*/
		template<class T> vector<T> apply_to_each_term(vector<T> const& v,
			UnaryOperation<T> const& un_op);

		template<class T> vector<T> apply_to_each_term(vector<T> const& v,
			UnaryOperation<T> const& un_op)
		{
			vector<T> res(v.size());
			for (std::size_t i = 0; i < v.size(); ++i)
				res[i] = un_op(v[i]);
			return res;
		}

		/*! \brief Selects a member of std::pair
		\param p A pair
		\param index 0 for first member, 1 for second, error otherwise
		\return Either first or second members of pair
		*/
		template<class T> T pair_member(const std::pair<T, T>& p, int index);

		template<class T> T pair_member(const std::pair<T, T>& p, int index)
		{
			assert((0 == index) || (1 == index));
			return 0 == index ? p.first : p.second;
		}

		/*! \brief Sets a member of std::pair
		\param p Pair
		\param index 0 for first, 1 for second, error otherwise
		\param val Value to be written
		*/
		template<class T> void set_pair_member(std::pair<T, T>& p, int index, const T& val);

		template<class T> void set_pair_member(std::pair<T, T>& p, int index, const T& val)
		{
			assert((0 == index) || (1 == index));
			if (0 == index)
				p.first = val;
			else
				p.second = val;
		}

		/*! \brief Performs type casting for an entire vector
		\param source Source vector
		\return Source vector converted to new type
		*/
		template<class T, class S> vector<T> list_static_cast(const vector<S>& source);

		template<class T, class S> vector<T> list_static_cast(const vector<S>& source)
		{
			vector<T> res(source.size());
			for (std::size_t i = 0; i < res.size(); ++i)
				res[i] = static_cast<T>(source[i]);
			return res;
		}

		/*! \brief Inserts all elements from one vector to the end of another
		\param subject Vector that will be modifies
		\param addendum Vector that will be added
		*/
		template<class T> void insert_all_to_back(vector<T>& subject, const vector<T>& addendum);

		template<class T> void insert_all_to_back(vector<T>& subject, const vector<T>& addendum)
		{
			if (!addendum.empty())
				subject.insert(subject.end(), addendum.begin(), addendum.end());
		}
		/*!
		\brief Binary search that return iterator of found object or end if not found
		\param begin begin iterator
		\param end end iterator
		\param val The value to search for
		\return iterator of found object or end if not found
		*/
		template<class Iter, class T>
		Iter binary_find(Iter begin, Iter end, T val);

		template<class Iter, class T>
		Iter binary_find(Iter begin, Iter end, T val)
		{
			// Finds the lower bound in at most log(last - first) + 1 comparisons
			Iter i = std::lower_bound(begin, end, val);

			if (i != end && !(val < *i))
				return i; // found
			else
				return end; // not found
		}

		template<typename T> size_t binary_index_find(std::vector<T> const& data, T const& key)
		{
			auto it = std::lower_bound(data.begin(), data.end(), key);
			if(it == data.end())
				throw MadVoro::Exception::MadVoroException("Key not found in binary_index_find");
			size_t const index = static_cast<size_t>(it - data.begin());
			if(data[index] != key)
				throw MadVoro::Exception::MadVoroException("Key not equal found in binary_index_find");
			return index;
		}

		/*! \brief Checks for existence and retrieves entry from flat map
		\param data Data vector
		\param key Key to look for
		\param keys Sorted vector of keys
		\return Value corresponding to key
		*/
		template<class S, class T> typename vector<T>::const_reference safe_retrieve
		(vector<T> const& data, vector<S> const& keys, const S& key);

		template<class S, class T> typename vector<T>::const_reference safe_retrieve
		(vector<T> const& data, vector<S> const& keys, const S& key)
		{
			assert(data.size() == keys.size());
			typename vector<S>::const_iterator it = binary_find(keys.begin(), keys.end(), key);
			assert(it != keys.end());
			std::size_t index = static_cast<std::size_t>(it - keys.begin());
			return data.at(index);
		}

		/*! \brief Retrieve a value from the container
		\param data_begin First iterator
		\param key_begin First key
		\param key_end Last key
		\param key Desired key
		\return Iterator to value
		*/
		template<typename It,typename It2, class S>
		const It safe_retrieve(const It data_begin, const It2 key_begin, const It2 key_end, S const& key)
		{
			const It2 index = binary_find(key_begin, key_end, key);
			assert(index != key_end);
			return data_begin + (index - key_begin);
		}


		/*! \brief Checks for existence and retrieves entry from flat map
		\param data Data vector
		\param key Key to look for
		\param keys Sorted vector of keys
		\return Value corresponding to key
		*/
		template<class S, class T> typename vector<T>::reference safe_retrieve
		(vector<T> &data, vector<S> const& keys, const S& key);

		template<class S, class T> typename vector<T>::reference safe_retrieve
		(vector<T> &data, vector<S> const& keys, const S& key)
		{
			assert(data.size() == keys.size());
			typename vector<S>::const_iterator it = binary_find(keys.begin(), keys.end(), key);
			assert(it != keys.end());
			return data[static_cast<std::size_t>(it - keys.begin())];
		}

		/*! \brief Checks for existence and retrieves entry from flat map
		\param m flat map
		\param s Key
		\return Value corresponding to key
		*/
		template<class S, class T> const T& safe_retrieve
		(const boost::container::flat_map<S, T>& m,
			const S& s);

		template<class S, class T> const T& safe_retrieve
		(const boost::container::flat_map<S, T>& m,
			const S& s)
		{
			typename boost::container::flat_map<S, T>::const_iterator it =
				m.find(s);
			assert(it != m.end());
			return it->second;
		}

		/*! \brief Non constant version of safe retrieve
		\param m flat map
		\param s Key
		\return Value corresponding to key
		*/
		template<class S, class T> T& safe_retrieve
		(boost::container::flat_map<S, T>& m,
			const S& s);

		template<class S, class T> T& safe_retrieve
		(boost::container::flat_map<S, T>& m,
			const S& s)
		{
			typename boost::container::flat_map<S, T>::iterator it =
				m.find(s);
			assert(it != m.end());
			return it->second;
		}

		/*! \brief Fast approximate sqrt
		\param x The value to calculate the sqrt of
		\return Sqrt(x)
		*/
		#ifdef __INTEL_COMPILER
		#pragma omp declare simd
		#endif
		double fastsqrt(double x);

		/**
		 * Function to perform 2D interpolation on a given table.
		 *
		 * @param x The x-coordinate for interpolation.
		 * @param y The y-coordinate for interpolation.
		 * @param x_vec The vector of x-coordinates for the table.
		 * @param y_vec The vector of y-coordinates for the table.
		 * @param data The 2D data table.
		 * @param x_vec_high_slope The slope for x-values greater than the maximum x-value in the table.
		 * @param slope_length The number of points to use for calculating the slope.
		 *
		 * @return The interpolated value at the given (x, y) coordinates.
		 *
		 * @throws MadVoroException If the size of x_vec or y_vec is less than or equal to slope_length.
		 */
		double Interpolate2DTable(double const x, double const y, std::vector<double> const& x_vec, std::vector<double> const& y_vec, std::vector<std::vector<double>> const& data,
				double const x_vec_high_slope = 0, size_t const slope_length = 5);
	}
}

#endif // UTILS_HPP
