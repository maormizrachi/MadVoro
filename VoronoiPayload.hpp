#ifndef VORONOI_PAYLOAD_HPP
#define VORONOI_PAYLOAD_HPP

#ifdef MADVORO_WITH_MPI
#include <mpi_utils/serialize/Serializable.hpp>
#include <mpi_utils/serialize/Serializer.hpp>
#endif

namespace MadVoro {

template <typename PointT>
struct VoronoiPayload
#ifdef MADVORO_WITH_MPI
    : public Serializable
#endif
{
    using coord_type = typename PointT::coord_type;

    coord_type radius;
    PointT CM;

    VoronoiPayload() : radius(0), CM() {}
    VoronoiPayload(coord_type r, const PointT &cm) : radius(r), CM(cm) {}

#ifdef MADVORO_WITH_MPI
    size_t dump(Serializer *serializer) const override
    {
        size_t count = serializer->insert(radius);
        count += CM.dump(serializer);
        return count;
    }

    size_t load(const Serializer *serializer, size_t byteOffset) override
    {
        size_t bytes = 0;
        bytes += serializer->extract(radius, byteOffset);
        bytes += CM.load(serializer, byteOffset + bytes);
        return bytes;
    }
#endif // MADVORO_WITH_MPI
};

} // namespace MadVoro

#endif // VORONOI_PAYLOAD_HPP
