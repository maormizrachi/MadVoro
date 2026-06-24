#ifndef DELAUNAY3D_HPP
#define DELAUNAY3D_HPP 1

#include "3D/elementary/Vector3D.hpp"
#include "3D/elementary/Tetrahedron.hpp"
#include <string>
#include <vector>
#include <set>
#include <array>
#include <boost/container/flat_set.hpp>

using std::vector;
using std::string;
//using std::stack;

//! \brief A three dimensional delauny triangulation
class Delaunay3D
{
public:
  //! \brief List of tetrahedra
  vector<Tetrahedron> tetras_;
  //! \brief List of mesh generating points
  vector<Vector3D> points_;
  //! \brief List of empty tetrahedra
  boost::container::flat_set<size_t> empty_tetras_;
  //! \brief Length of list of real points (the first `Norg_` points from `points_` are actually mine, the rest are ghosts)
  std::size_t Norg_; 
  //! \brief List of outside neighbours
  std::size_t outside_neighbor_;
  std::vector<std::pair<size_t, Tetrahedron>> changed_tetras_;
  std::vector<size_t> newTetras_;

  Delaunay3D();

  /*! \brief Copy constructor
    \param other Source
   */
  Delaunay3D(Delaunay3D const& other);

  /*! \brief Copy assingment
    \param other Source
    \return Reference to new object
   */
  Delaunay3D& operator=(Delaunay3D const& other);

  ~Delaunay3D();

  /*! \brief Build triangulation
    \param points Mesh generating points
    \param maxv Bounding point
    \param minv Bounding point
    \param order Sorted list of indices
   */ 
  void Build(vector<Vector3D> const& points,Vector3D const& maxv,Vector3D const& minv,std::vector<size_t> &order);

  /*! \brief Add new points and rebuild
    \param points New points
   */
  void BuildExtra(vector<Vector3D> const& points);

  /*! \brief Dump data to file
    \param filename Name of output file
   */
  void output(string const& filename) const;

  /*! \brief Validate input
    \return True if input is valid
   */
  bool CheckCorrect(void);

  //! \brief Clear data
  void Clean(void);
private:
  void InsertPoint(std::size_t index);
  std::size_t Walk(std::size_t point, std::size_t first_guess);
  void flip14(std::size_t point,std::size_t tetra);
  void flip23(std::size_t tetra0, std::size_t tetra1,std::size_t location0,bool flat_check);
  void flip32(std::size_t tetra0, std::size_t tetra1, std::size_t location0, std::size_t shared_loction, bool flat_check);
  void flip44(std::size_t tetra0, std::size_t tetra1, std::size_t location0, std::size_t neigh0, std::size_t neigh1);
  void FindFlip(std::size_t tetra0, std::size_t tetra1, std::size_t p, size_t p_loc, size_t other_point_loc);
  void ExactFlip(std::size_t tetra0, std::size_t tetra1, std::size_t p);
  std::size_t FindThirdNeighbor(std::size_t tetra0, std::size_t tetra1);

  std::array<Vector3D, 3> b3_temp_,b3_temp2_;
  std::array<Vector3D, 4> b4_temp_;
  std::array<Vector3D, 5> b5_temp_;
  //  std::array<std::size_t, 4> b4s_temp_, b4s_temp2_;
  std::array<std::size_t, 8> b8s_temp_;
  vector<std::size_t> to_check_;
  std::size_t last_checked_;	
  Tetrahedron tet_temp0_, tet_temp1_, newtet_;
};

#endif //DELAUNAY3D_HPP
