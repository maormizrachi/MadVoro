#ifndef MADVORO_FACE3D_HPP
#define MADVORO_FACE3D_HPP

#include <vector>
#include <numeric>
#include <cstddef>
#include <boost/container/small_vector.hpp>
#include "PointOps.hpp"

namespace MadVoro {

using namespace MadVoro::fallback;

template <typename PointT>
using point_vec_v_t = boost::container::small_vector<PointT, 10>;

template <typename PointT>
class Face3D
{
public:
  point_vec_v_t<PointT> vertices;
  std::pair<std::size_t, std::size_t> neighbors;

  Face3D() : vertices(), neighbors() {}

  Face3D(point_vec_v_t<PointT> const& vert, std::size_t neighbor1, std::size_t neighbor2)
    : vertices(vert), neighbors(neighbor1, neighbor2) {}

  Face3D(Face3D const& other) = default;
  Face3D& operator=(Face3D const& other) = default;
  ~Face3D() = default;

  double GetArea() const
  {
    const PointT& ref = vertices[0];
    return std::inner_product(vertices.begin()+1,
                              vertices.end()-1,
                              vertices.begin()+2,
                              0.0,
                              [](double x, double y){ return x + y; },
                              [&ref](const PointT& u, const PointT& v)
                              { return 0.5 * fastabs(CrossProduct(u - ref, v - ref)); });
  }
};

template <typename PointT>
PointT calc_centroid(const Face3D<PointT>& face)
{
  const PointT& ref = face.vertices[0];
  return std::inner_product(face.vertices.begin()+1,
                            face.vertices.end()-1,
                            face.vertices.begin()+2,
                            PointT(0, 0, 0),
                            [](const PointT& x, const PointT& y){ return x + y; },
                            [&ref](const PointT& u, const PointT& v)
                            {
                              const double area = 0.5 * fastabs(CrossProduct(u - ref, v - ref));
                              return area * (u + v + ref) / 3.0;
                            });
}

} // namespace MadVoro

#endif // MADVORO_FACE3D_HPP
