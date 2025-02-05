#ifndef FACE_HPP
#define FACE_HPP 1

#include <vector>
#include <numeric>
#include <boost/container/small_vector.hpp>
#include "Vector3D.hpp"

//! \brief Container for small collection of points

namespace MadVoro
{
  typedef boost::container::small_vector<Vector3D, 10> point_vec_v;

  //! \brief Interface between two cells
  class Face
  {
  public:

    //! \brief Points at the ends of the edge
    point_vec_v vertices;

    //! \brief Neighboring cells
    std::pair<std::size_t,std::size_t> neighbors;
       
      // TODO: get neighbors methods

    /*! \brief Class constructor
      \param vert Position of the vertices
      \param neighbor1 Index of first neighbor cell
      \param neighbor2 Index of second neighbor cell
    */
    Face(const point_vec_v &vert, std::size_t neighbor1, std::size_t neighbor2);

    /*! \brief Assignment operator
      \param other Source
      \return Reference to new object
    */
    Face& operator=(const Face& other);

    Face(void);

    ~Face(void);

    /*! \brief Copy constructor
      \param other Source Face3D
    */
    Face(Face const& other);
  };
}

#endif	// FACE_HPP
