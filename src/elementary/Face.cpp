#include "Face.hpp"
 
MadVoro::Face::Face(point_vec_v const& vert,std::size_t neighbor1,std::size_t neighbor2):
  vertices(vert),neighbors(neighbor1,neighbor2) {}

MadVoro::Face& MadVoro::Face::operator=(const Face& other)
{
  vertices = other.vertices;
  neighbors = other.neighbors;
  return *this;
}

MadVoro::Face::Face(void): vertices(), neighbors() {}

MadVoro::Face::~Face(void)
{}

MadVoro::Face::Face(Face const& other):
  vertices(other.vertices),
  neighbors(other.neighbors) {}