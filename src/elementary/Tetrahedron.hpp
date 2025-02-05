#ifndef TETRAHEDRON_HPP
#define TETRAHEDRON_HPP 1

#include <array>

//points are ordered such as that the fourth point is above the plane defined by points 0 1 2 in a couter clockwise fashion
// neighbors are the tetra opposite to the triangle starting with the index of the vertice

namespace MadVoro
{
  /*! \brief A class describing tetrahedrons
  */
  class Tetrahedron
  {
  public:
    //! \brief Indices of vertices
    std::size_t points[4];
    //! \brief Indices of neighbours
    std::size_t neighbors[4];
    bool newTetra;
    
    Tetrahedron();

    /*! \brief Copy constructor
      \param other Source
    */
    Tetrahedron(Tetrahedron const& other);

    ~Tetrahedron();

    /*! \brief Copy assignment
      \param other Source
      \return Reference to new object
    */
    Tetrahedron& operator=(Tetrahedron const& other);
  };
}

#endif //TETRAHEDRON_HPP
