#include "Face3D.hpp"

using std::inner_product;
using namespace MadVoro;
 
MadVoro::Face3D::Face3D(point_vec_v const& vert,std::size_t neighbor1,std::size_t neighbor2):
  vertices(vert),neighbors(neighbor1,neighbor2) {}

Face3D& MadVoro::Face3D::operator=(const Face3D& other)
{
  vertices = other.vertices;
  neighbors = other.neighbors;
  return *this;
}

MadVoro::Face3D::Face3D(void): vertices(), neighbors() {}

MadVoro::Face3D::~Face3D(void)
{}

MadVoro::Face3D::Face3D(Face3D const& other):
  vertices(other.vertices),
  neighbors(other.neighbors) {}

double MadVoro::Face3D::GetArea(void) const
{
  const Point3D& ref = vertices[0];
  return inner_product(vertices.begin()+1,
		       vertices.end()-1,
		       vertices.begin()+2,
		       0,
		       [](const double x, const double y)
		       {return x+y;},
		       [&ref](const Point3D& u, const Point3D& v)
		       {return 0.5*fastabs(CrossProduct(u-ref,v-ref));});
}

Point3D calc_centroid(const Face3D& face)
{
  const Point3D& ref = face.vertices[0];
  return inner_product(face.vertices.begin()+1,
		       face.vertices.end()-1,
		       face.vertices.begin()+2,
		       Point3D(0,0,0),
		       [](const Point3D& x, const Point3D& y)
		       {return x+y;},
		       [&ref](const Point3D& u, const Point3D& v)
		       {
			 const double area = 0.5*fastabs
			   (CrossProduct(u-ref, v-ref));
			 return area*(u+v+ref)/3;
		       });
}
