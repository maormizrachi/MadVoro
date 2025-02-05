/*! \file Mat33.hpp
\brief A very simple class for a 3x3 matrix
\author Elad Steinberg
*/
#ifndef MAT3_HPP
#define MAT3_HPP 1

#include "Point3D.hpp"
#include "mpi/serialize/Serializer.hpp"

#define EPSILON 1e-12

namespace MadVoro
{

	//! \brief A 3x3 matrix
	template <typename T>
	class Mat33 
				#ifdef MADVORO_WITH_MPI
					: public MadVoro::MPI::Serializable
				#endif // MADVORO_WITH_MPI
	{
	private:
		T _data[3][3];

	public:
		//! \brief Return the element at (row, col)
	//! \param row Row
	//! \param col Column
	//! \return Entry
		inline T& at(int row, int col)
		{
			return _data[row][col];
		}

		//! \brief Return the element at (row, col)
	//! \param row Row
	//! \param col Column
	//! \return Entry
		inline const T& at(int row, int col) const
		{		
			return _data[row][col];
		}

		inline void Set(T d00, T d01, T d02, T d10, T d11, T d12, T d20, T d21, T d22)
		{
			_data[0][0] = d00;
			_data[0][1] = d01;
			_data[0][2] = d02;
			_data[1][0] = d10;
			_data[1][1] = d11;
			_data[1][2] = d12;
			_data[2][0] = d20;
			_data[2][1] = d21;
			_data[2][2] = d22;
		}

		inline void SetAt(T value, int i, int j)
		{
			_data[i][j] = value;
		}

		inline void AddAt(T value, int i, int j)
		{
			_data[i][j]+=value;
		}


		//! \brief Return the element at (row, col)
	//! \param row Row
	//! \param col Column
	//! \return Entry
		inline T& operator()(int row, int col) { return at(row, col); }

		//! \brief Return the element at (row, col)
	//! \param row Row
	//! \param col Column
	//! \return Entry
		inline const T& operator()(int row, int col) const { return at(row, col); }

		//! \brief Constructs a zeroed out matrix
		Mat33();

		//! \brief Constructs a matrix and fills all its values. An intializer_list is better, but C++ 11 isn't always supported
	/*! \param d00 Term in position 0,0
		\param d01 Term in position 0,1
		\param d02 Term in position 0,2
		\param d10 Term in position 1,0
		\param d11 Term in position 1,1
		\param d12 Term in position 1,2
		\param d20 Term in position 2,0
		\param d21 Term in position 2,1
		\param d22 Term in position 2,2
	*/
		Mat33(T d00, T d01, T d02,
			T d10, T d11, T d12,
			T d20, T d21, T d22);

		//! \brief Returns the matrix's determinant
	//! \return Determinant of the matrix
		T determinant() const;

	//! \brief Return the inverse matrix
	//! \return Inverse matrix
		Mat33<T> inverse()const;

	//! \brief Returns the transpose matrix
	//! \return Transposed matrix
		Mat33<T> transpose()const;

		/*! \brief Addition
		\param A Matrix to be added
		\return Reference to sum
		*/
	#ifdef __INTEL_COMPILER
	#pragma omp declare simd
	#endif
		Mat33& operator+=(Mat33 const& A);

		/*! \brief Subtraction
		\param A Matrix to be subtracted
		\return Difference
		*/
	#ifdef __INTEL_COMPILER
	#pragma omp declare simd
	#endif
		Mat33& operator-=(Mat33 const& A);

		/*! \brief Assignment operator
		\param A Matrix to be copied
		\return The assigned value
		*/
	#ifdef __INTEL_COMPILER
	#pragma omp declare simd
	#endif
		Mat33& operator=(Mat33 const& A);
		
		/*! \brief Scalar product
		\param s Scalar
		\return Reference to the Matrix multiplied by scalar
		*/
	#ifdef __INTEL_COMPILER
	#pragma omp declare simd
	#endif
		Mat33& operator*=(double s);

		/*! \brief Compare 3D-Vectors (up to an arbitrary precision)
		\param A Matrix to be compared to
		\return True/False - according to the comparison results.
		*/
		bool operator==(Mat33 const& A) const;


		/*! \brief Computes the second invariant of the matrix (J2)
		\return double - J2
		*/
		T J2() const;


	#ifdef __INTEL_COMPILER
	#pragma omp declare simd
	#endif
		~Mat33(void) {}


	#ifdef MADVORO_WITH_MPI
		force_inline size_t dump(MadVoro::MPI::Serializer *serializer) const override
		{
			size_t bytes = 0;
			bytes += serializer->insert(this->_data[0][0]);
			bytes += serializer->insert(this->_data[0][1]);
			bytes += serializer->insert(this->_data[0][2]);
			bytes += serializer->insert(this->_data[1][0]);
			bytes += serializer->insert(this->_data[1][1]);
			bytes += serializer->insert(this->_data[1][2]);
			bytes += serializer->insert(this->_data[2][0]);
			bytes += serializer->insert(this->_data[2][1]);
			bytes += serializer->insert(this->_data[2][2]);
			return bytes;
		}

		force_inline size_t load(const MadVoro::MPI::Serializer *serializer, std::size_t byteOffset) override
		{
			size_t bytes = 0;
			bytes += serializer->extract(this->_data[0][0], byteOffset);
			bytes += serializer->extract(this->_data[0][1], byteOffset + bytes);
			bytes += serializer->extract(this->_data[0][2], byteOffset + bytes);
			bytes += serializer->extract(this->_data[1][0], byteOffset + bytes);
			bytes += serializer->extract(this->_data[1][1], byteOffset + bytes);
			bytes += serializer->extract(this->_data[1][2], byteOffset + bytes);
			bytes += serializer->extract(this->_data[2][0], byteOffset + bytes);
			bytes += serializer->extract(this->_data[2][1], byteOffset + bytes);
			bytes += serializer->extract(this->_data[2][2], byteOffset + bytes);
			return bytes;
		}
	#endif // MADVORO_WITH_MPI
	};

	/*! \brief Term by term addition
		\param A1 First Matrix
		\param A2 Second Matrix
		\return Sum
		*/
		template <typename T>
		Mat33<T> operator+(Mat33<T> const& A1, Mat33<T> const& A2);

		/*! \brief Term by term subtraction
		\param A1 First Matrix
		\param A2 Second Matrix
		\return Difference
		*/
		template <typename T>
		Mat33<T> operator-(Mat33<T> const& A1, Mat33<T> const& A2);

		/*! \brief Scalar product
		\param d Scalar
		\param A Matrix
		\return Three dimensional Matrix
		*/
		template <typename T>
		Mat33<T> operator*(double d, Mat33<T> const& A);

		/*! \brief Matrix Matrix product
		\param A1 Matrix
		\param A2 Matrix
		\return Three dimensional Matrix
		*/
		template <typename T>
		Mat33<T> operator*(Mat33<T> const& A1, Mat33<T> const& A2);

		/*! \brief Matrix Vector product
		\param A Matrix
		\param v Vector
		\return Three dimensional Matrix
		*/
		template <typename T>
		Point3D operator*(Mat33<T> const& A,  Point3D const& v);

		/*! \brief Scalar product
		\param A Matrix
		\param d Scalar
		\return Three dimensional Matrix
		*/
		template <typename T>
		Mat33<T> operator*(Mat33<T> const& A, double d);

		/*! \brief Scalar division
		\param A Matrix
		\param d Scalar
		\return Three dimensional Matrix
		*/
		template <typename T>
		Mat33<T> operator/(Mat33<T> const& A, double d);

		/*! \brief Sum_{i,j} A1ij*A2ij
		\param A1 Matrix
		\param A2 Matrix
		\return double
		*/
		template <typename T>
		T operator%(Mat33<T> const& A1, Mat33<T> const& A2);


		/*! \brief retrun the deviatoric part
		\param A Matrix
		\return A-1/3tr(A)*I
		*/
		template <typename T>
		Mat33<T> deviator(Mat33<T>const & A);

	template <typename T>
	inline T Mat33<T>::determinant() const
	{
		return at(0, 0)*(at(1, 1)*at(2, 2) - at(1, 2)*at(2, 1)) + at(0, 1)*(at(1, 2)*at(2, 0) - at(1, 0)*at(2, 2))
			+ at(0, 2)*(at(1, 0)*at(2, 1) - at(1, 1)*at(2, 0));
	}

	template<typename T>
	inline Mat33<T> Mat33<T>::inverse() const
	{
		Mat33<T> temp;
		double det = determinant();
		double detInverse = 1 / det;
		temp._data[0][0] = (_data[1][1] * _data[2][2] - _data[2][1] * _data[1][2]) * detInverse;
		temp._data[0][1] = (_data[0][2] * _data[2][1] - _data[0][1] * _data[2][2]) * detInverse;
		temp._data[0][2] = (_data[0][1] * _data[1][2] - _data[0][2] * _data[1][1]) * detInverse;
		temp._data[1][0] = (_data[1][2] * _data[2][0] - _data[1][0] * _data[2][2]) * detInverse;
		temp._data[1][1] = (_data[0][0] * _data[2][2] - _data[0][2] * _data[2][0]) * detInverse;
		temp._data[1][2] = (_data[0][2] * _data[1][0] - _data[0][0] * _data[1][2]) * detInverse;
		temp._data[2][0] = (_data[1][0] * _data[2][1] - _data[1][1] * _data[2][0]) * detInverse;
		temp._data[2][1] = (_data[0][1] * _data[2][0] - _data[0][0] * _data[2][1]) * detInverse;
		temp._data[2][2] = (_data[1][1] * _data[0][0] - _data[0][1] * _data[1][0]) * detInverse;
		return temp;
	}

	template<typename T>
	inline Mat33<T> Mat33<T>::transpose() const
	{
		Mat33<T> res;
		for (size_t i = 0; i < 3; ++i)
			for (size_t j = 0; j < 3; ++j)
				res._data[i][j] = _data[j][i];
		return res;
	}

	template <typename T>
	inline Mat33<T>::Mat33(T d00, T d01, T d02, 
		T d10, T d11, T d12, 
		T d20, T d21, T d22)
	{
		_data[0][0] = d00;
		_data[0][1] = d01;
		_data[0][2] = d02;
		_data[1][0] = d10;
		_data[1][1] = d11;
		_data[1][2] = d12;
		_data[2][0] = d20;
		_data[2][1] = d21;
		_data[2][2] = d22;
	}

	template<typename T>
	Mat33<T>::Mat33()
	{
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				_data[i][j] = 0;
	}

	template <typename T>
	inline T Mat33<T>::J2() const
	{
		return _data[0][0]*_data[0][0] + _data[0][1]*_data[0][1] + _data[0][2]*_data[0][2] 
			+ _data[1][0]*_data[1][0] + _data[1][1]*_data[1][1] + _data[1][2]*_data[1][2]
			+ _data[2][0]*_data[2][0] + _data[2][1]*_data[2][1] + _data[2][2]*_data[2][2];
	}

	template <typename T>
	Mat33<T> operator+(Mat33<T> const& A1, Mat33<T> const& A2)
	{
		Mat33<T> res(A1.at(0,0)+A2.at(0,0), A1.at(0,1)+A2.at(0,1), A1.at(0,2)+A2.at(0,2),
					A1.at(1,0)+A2.at(1,0), A1.at(1,1)+A2.at(1,1), A1.at(1,2)+A2.at(1,2),
					A1.at(2,0)+A2.at(2,0), A1.at(2,1)+A2.at(2,1), A1.at(2,2)+A2.at(2,2));
		return res;
	}

	template <typename T>
	Mat33<T> operator-(Mat33<T> const& A1, Mat33<T> const& A2)
	{
		Mat33<T> res(A1.at(0,0)-A2.at(0,0), A1.at(0,1)-A2.at(0,1), A1.at(0,2)-A2.at(0,2),
					A1.at(1,0)-A2.at(1,0), A1.at(1,1)-A2.at(1,1), A1.at(1,2)-A2.at(1,2),
					A1.at(2,0)-A2.at(2,0), A1.at(2,1)-A2.at(2,1), A1.at(2,2)-A2.at(2,2));
		return res;
	}

	template <typename T>
	Mat33<T> operator*(double d, Mat33<T> const& A)
	{
		Mat33<T> res(A.at(0,0)*d, A.at(0,1)*d, A.at(0,2)*d,
					A.at(1,0)*d, A.at(1,1)*d, A.at(1,2)*d,
					A.at(2,0)*d, A.at(2,1)*d, A.at(2,2)*d);
		return res;
	}

	template <typename T>
	Mat33<T> operator*(Mat33<T> const& A, double d)
	{
		Mat33<T> res(A.at(0,0)*d, A.at(0,1)*d, A.at(0,2)*d,
					A.at(1,0)*d, A.at(1,1)*d, A.at(1,2)*d,
					A.at(2,0)*d, A.at(2,1)*d, A.at(2,2)*d);
		return res;
	}

	template <typename T>
	Mat33<T> operator/(Mat33<T> const& A, double d)
	{
		Mat33<T> res(A.at(0,0)/d, A.at(0,1)/d, A.at(0,2)/d,
					A.at(1,0)/d, A.at(1,1)/d, A.at(1,2)/d,
					A.at(2,0)/d, A.at(2,1)/d, A.at(2,2)/d);
		return res;
	}

	template <typename T>
	T operator%(Mat33<T> const &A1, Mat33<T> const &A2)
	{
		return A1(0,0)*A2(0,0)+A1(0,1)*A2(0,1)+A1(0,2)*A2(0,2)+A1(1,0)*A2(1,0)+A1(1,1)*A2(1,1)+A1(1,2)*A2(1,2)+A1(2,0)*A2(2,0)+A1(2,1)*A2(2,1)+A1(2,2)*A2(2,2);
	}

	#ifdef __INTEL_COMPILER
	#pragma omp declare simd
	#endif
	template <typename T>
	Mat33<T>& Mat33<T>::operator+=(Mat33<T> const& A)
	{
		for (int i=0; i<3; ++i)
		{
			for (int j=0; j<3; ++j)
			{
				_data[i][j] += A(i, j);
			}
		}
		return *this;
	}

	#ifdef __INTEL_COMPILER
	#pragma omp declare simd
	#endif
	template <typename T>
	Mat33<T>& Mat33<T>::operator-=(Mat33<T> const& A)
	{
		for (int i=0; i<3; ++i)
		{
			for (int j=0; j<3; ++j)
			{
				_data[i][j] -= A(i, j);
			}
		}
		return *this;
	}

	#ifdef __INTEL_COMPILER
	#pragma omp declare simd
	#endif
	template <typename T>
	Mat33<T>& Mat33<T>::operator*=(double d)
	{
		for (int i=0; i<3; ++i)
		{
			for (int j=0; j<3; ++j)
			{
				_data[i][j] *= d;
			}
		}
		return *this;
	}

	#ifdef __INTEL_COMPILER
	#pragma omp declare simd
	#endif
	template <typename T>
	Mat33<T>& Mat33<T>::operator=(Mat33<T> const& A)
	{
		for (int i =0; i<3; ++i)
		{
			for (int j=0; j<3; ++j)
			{
				_data[i][j] = A.at(i, j);
			}
		}
		return *this;
	}

	template <typename T>
	bool Mat33<T>::operator==(Mat33<T> const& A) const
	{
		bool res;
		for (int i =0; i<3; ++i)
		{
			for (int j=0; j<3; ++j)
			{
				res = res && std::abs(_data[i][j] - A.at(i, j)) < EPSILON;
			}
		}
		return res;
	}

	template <typename T>
	Point3D operator*(Mat33<T> const& A, Point3D const& v)
	{
		Point3D res(A.at(0,0)*v[0]+A.at(0, 1)*v[1]+A.at(0,2)*v[2],
					A.at(1,0)*v[0]+A.at(1, 1)*v[1]+A.at(1,2)*v[2],
					A.at(2,0)*v[0]+A.at(2, 1)*v[1]+A.at(2,2)*v[2]);

		return res;
	}


	template <typename T>
	Mat33<T> operator*(Mat33<T> const& A1, Mat33<T> const& A2)
	{
		Mat33<T> res(A1(0, 0)*A2(0, 0) + A1(0, 1)*A2(1, 0) + A1(0, 2)*A2(2, 0),	   A1(0, 0)*A2(0, 1) + A1(0, 1)*A2(1, 1) + A1(0, 2)*A2(2, 1),      A1(0, 0)*A2(0, 2) + A1(0, 1)*A2(1, 2) + A1(0, 2)*A2(2, 2),
					A1(1, 0)*A2(0, 0) + A1(1, 1)*A2(1, 0) + A1(1, 2)*A2(2, 0),    A1(1, 0)*A2(0, 1) + A1(1, 1)*A2(1, 1) + A1(1, 2)*A2(2, 1),      A1(1, 0)*A2(0, 2) + A1(1, 1)*A2(1, 2) + A1(1, 2)*A2(2, 2),
					A1(2, 0)*A2(0, 0) + A1(2, 1)*A2(1, 0) + A1(2, 2)*A2(2, 0),    A1(2, 0)*A2(0, 1) + A1(2, 1)*A2(1, 1) + A1(2, 2)*A2(2, 1),      A1(2, 0)*A2(0, 2) + A1(2, 1)*A2(1, 2) + A1(2, 2)*A2(2, 2));
		return res;
	}

	template <typename T>
	Mat33<T> deviator(Mat33<T> const& A)
	{
		T trace3 = (A(0,0)+A(1,1)+A(2,2))*1/3.;
		Mat33<T> res(A(0,0)-trace3, A(0,1),A(0,2), A(1,0), A(1,1)-trace3, A(1,2), A(2,0), A(2,1), A(2,2)-trace3);
		return res;
	}
}

#endif // MAT33_HPP
