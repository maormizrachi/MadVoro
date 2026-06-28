#ifndef RANGE_QUERY_DATA
#define RANGE_QUERY_DATA

#include <cstddef>

#ifdef MADVORO_WITH_MPI
    #include <mpi_utils/serialize/Serializer.hpp>
#endif // MADVORO_WITH_MPI

template <typename PointT>
struct RangeQueryData 
#ifdef MADVORO_WITH_MPI
    : public Serializable
#endif // MADVORO_WITH_MPI
{
    using coord_type = typename PointT::coord_type;

    size_t pointIdx;
    PointT center;
    coord_type radius;

    RangeQueryData(size_t pointIdx, const PointT &center, coord_type radius):
        pointIdx(pointIdx), center(center), radius(radius)
    {};

    RangeQueryData(): pointIdx(0), center(PointT()), radius(0){};
    
    #ifdef MADVORO_WITH_MPI
        inline size_t dump(Serializer *serializer) const override
        {
            size_t bytes = 0;
            bytes += serializer->insert(this->pointIdx);
            bytes += serializer->insert(this->center);
            bytes += serializer->insert(this->radius);
            return bytes;
        }

        inline size_t load(const Serializer *serializer, size_t byteOffset) override
        {
            size_t bytes = 0;
            bytes += serializer->extract(this->pointIdx, byteOffset);
            bytes += serializer->extract(this->center, byteOffset + bytes);
            bytes += serializer->extract(this->radius, byteOffset + bytes);
            return bytes;
        }
    #endif // MADVORO_WITH_MPI
};

#endif // RANGE_QUERY_DATA