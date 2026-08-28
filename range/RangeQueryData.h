#ifndef RANGE_QUERY_DATA
#define RANGE_QUERY_DATA

#include <cstddef>
#include <cmath>
#include <limits>

#ifdef MADVORO_WITH_MPI
    #include <mpi_utils/serialize/Serializer.hpp>
#endif // MADVORO_WITH_MPI

static constexpr int NUM_IMAGE_CODES = 27;
static constexpr int ZERO_IMAGE_CODE = 13;

template <typename PointT>
inline int ImageCode(const PointT &translation)
{
    auto signComponent = [](typename PointT::coord_type v) -> int
    {
        if(v > typename PointT::coord_type(0))
        {
            return 1;
        }
        if(v < typename PointT::coord_type(0))
        {
            return -1;
        }
        return 0;
    };
    int sx = signComponent(translation.x);
    int sy = signComponent(translation.y);
    int sz = signComponent(translation.z);
    return (sx + 1) * 9 + (sy + 1) * 3 + (sz + 1);
}

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
    PointT imageTranslation;

    RangeQueryData(size_t pointIdx, const PointT &center, coord_type radius):
        pointIdx(pointIdx), center(center), radius(radius), imageTranslation(PointT())
    {};

    RangeQueryData(): pointIdx(0), center(PointT()), radius(0), imageTranslation(PointT()){};
    
    #ifdef MADVORO_WITH_MPI
        inline size_t dump(Serializer *serializer) const override
        {
            size_t bytes = 0;
            bytes += serializer->insert(this->pointIdx);
            bytes += serializer->insert(this->center);
            bytes += serializer->insert(this->radius);
            bytes += serializer->insert(this->imageTranslation);
            return bytes;
        }

        inline size_t load(const Serializer *serializer, size_t byteOffset) override
        {
            size_t bytes = 0;
            bytes += serializer->extract(this->pointIdx, byteOffset);
            bytes += serializer->extract(this->center, byteOffset + bytes);
            bytes += serializer->extract(this->radius, byteOffset + bytes);
            bytes += serializer->extract(this->imageTranslation, byteOffset + bytes);
            return bytes;
        }
    #endif // MADVORO_WITH_MPI
};

#endif // RANGE_QUERY_DATA