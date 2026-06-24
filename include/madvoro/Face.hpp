#ifndef MADVORO_FACE_HPP
#define MADVORO_FACE_HPP

#include <cstddef>
#include <utility>
#include <initializer_list>
#include <boost/container/small_vector.hpp>

namespace MadVoro
{
  template<typename Vec3>
  class Face
  {
  public:
    using point_vec_v = boost::container::small_vector<Vec3, 10>;

    point_vec_v vertices;
    std::pair<std::size_t, std::size_t> neighbors;

    Face() : vertices(), neighbors() {}

    Face(const point_vec_v &vert, std::size_t neighbor1, std::size_t neighbor2)
      : vertices(vert), neighbors(neighbor1, neighbor2) {}

    Face(std::initializer_list<Vec3> vert)
      : vertices(vert.begin(), vert.end()), neighbors(0, 0) {}

    Face(const point_vec_v &vert)
      : vertices(vert), neighbors(0, 0) {}

    Face(const Face &other) = default;
    Face& operator=(const Face &other) = default;
    ~Face() = default;
  };
}

#endif // MADVORO_FACE_HPP
