#ifndef FACE3D_HPP
#define FACE3D_HPP 1

#include <vector>
#include <numeric>
#include <boost/container/small_vector.hpp>
#include "Point3D.hpp"
//! \brief Container for small collection of points

namespace MadVoro
{
  //! \brief Interface between two cells
  class Face3D
  {
  public:

    typedef boost::container::small_vector<Point3D, 10> point_vec_v;

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
    Face3D(point_vec_v const& vert,std::size_t neighbor1,std::size_t neighbor2);

    /*! \brief Assignment operator
      \param other Source
      \return Reference to new object
    */
    Face3D& operator=(const Face3D& other);

    Face3D(void);

    ~Face3D(void);

    /*! \brief Copy constructor
      \param other Source Face3D
    */
    Face3D(Face3D const& other);

    /*! \brief Returns the area of the face
      \return Length
    */
    double GetArea(void) const;
  };

  /*! \brief Calculates the centroid of aa face
    \param face The face
    \return Centroid
  */
  Point3D calc_centroid(const Face3D& face);
}

#endif	// FACE3D_HPP
