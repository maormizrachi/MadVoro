#ifndef MADVORO_TYPES_HPP
#define MADVORO_TYPES_HPP

#include <cstddef>
#include <boost/container/small_vector.hpp>
#include <boost/container/flat_map.hpp>

namespace MadVoro {

typedef boost::container::small_vector<std::size_t, 24> face_vec;
typedef boost::container::small_vector<std::size_t, 8> point_vec;
typedef boost::container::small_vector<std::size_t, 40> tetra_vec;
using AllPointsMap = boost::container::flat_map<std::size_t, std::size_t>;

} // namespace MadVoro

#endif // MADVORO_TYPES_HPP
