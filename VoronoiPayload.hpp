#ifndef VORONOI_PAYLOAD_HPP
#define VORONOI_PAYLOAD_HPP

#include "3D/elementary/Vector3D.hpp"
#include <mpi_utils/serialize/Serializable.hpp>
#include <mpi_utils/serialize/Serializer.hpp>

struct VoronoiPayload : public Serializable
{
    double radius;
    Vector3D CM;

    VoronoiPayload() : radius(0), CM() {}
    VoronoiPayload(double r, const Vector3D &cm) : radius(r), CM(cm) {}

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
};

#endif // VORONOI_PAYLOAD_HPP
