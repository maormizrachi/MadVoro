#ifndef MADVORO_GEOMETRY_KERNELS_MIRROR_CM_HPP
#define MADVORO_GEOMETRY_KERNELS_MIRROR_CM_HPP

#include "../../elementary/PointOps.hpp"

namespace MadVoro
{
namespace Geometry
{

//! \brief Mirror image of the real cell's centroid across a rigid boundary face.
template<typename PointT>
inline PointT MirroredCM(const PointT &first, const PointT &second, const PointT &realCM, const PointT &faceVertex)
{
    using namespace MadVoro::fallback;
    const PointT normal = normalize(first - second);
    return realCM - 2 * normal * ScalarProd(normal, realCM - faceVertex);
}

} // namespace Geometry
} // namespace MadVoro

#endif // MADVORO_GEOMETRY_KERNELS_MIRROR_CM_HPP
