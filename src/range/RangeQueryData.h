#ifndef RANGE_QUERY_DATA
#define RANGE_QUERY_DATA

#include "utils/point/3DPoint.hpp"

namespace MadVoro
{
    struct RangeQueryData 
                        #ifdef MADVORO_WITH_MPI
                            : public MadVoro::MPI::Serializable
                        #endif // MADVORO_WITH_MPI
    {
        size_t pointIdx;
        _3DPoint center;
        typename _3DPoint::coord_type radius;

        RangeQueryData(size_t pointIdx, const _3DPoint &center, typename _3DPoint::coord_type radius):
            pointIdx(pointIdx), center(center), radius(radius)
        {};

        RangeQueryData(): pointIdx(0), center(_3DPoint()), radius(0){};
        
        #ifdef MADVORO_WITH_MPI
            inline size_t dump(MadVoro::MPI::Serializer *serializer) const override
            {
                size_t bytes = 0;
                bytes += serializer->insert(this->pointIdx);
                bytes += serializer->insert(this->center);
                bytes += serializer->insert(this->radius);
                return bytes;
            }

            inline size_t load(const MadVoro::MPI::Serializer *serializer, size_t byteOffset) override
            {
                size_t bytes = 0;
                bytes += serializer->extract(this->pointIdx, byteOffset);
                bytes += serializer->extract(this->center, byteOffset + bytes);
                bytes += serializer->extract(this->radius, byteOffset + bytes);
                return bytes;
            }
        #endif // MADVORO_WITH_MPI
    };
}

#endif // RANGE_QUERY_DATA