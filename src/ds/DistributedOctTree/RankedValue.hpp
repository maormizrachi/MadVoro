#ifndef RANKED_VALUE_HPP
#define RANKED_VALUE_HPP

#include <boost/container/small_vector.hpp>

#define UNDEFINED_OWNER -1
#define PATH_END_DIRECTION (-1)

namespace MadVoro
{
    namespace DataStructure
    {
        template<typename T>
        using BoundingBox = MadVoro::Geometry::BoundingBox<T>;

        template<typename T, int max_ranks_per_leaf, int max_directions_size>
        struct OctTreeRankedValue
        {
            using coord_type = typename T::coord_type;
            using Raw_type = typename is_raw_type_defined<T>::type;

            T value;
            boost::container::small_vector<int, max_ranks_per_leaf> owners; // ranks of the owners (if several)
            direction_t directions[max_directions_size]; // where it is in the owner's tree (directions)

            OctTreeRankedValue(const T &value = T(), int owner = UNDEFINED_OWNER): value(value)
            {
                size_t numOwners = (owner == UNDEFINED_OWNER)? 0 : 1;
                this->owners.resize(numOwners);
                if(numOwners)
                {
                    this->owners[0] = owner;
                }
                this->directions[0] = PATH_END_DIRECTION;
            };

            OctTreeRankedValue(const OctTreeRankedValue &other)
            {
                (*this) = other;
            }

            OctTreeRankedValue &operator=(const OctTreeRankedValue &other)
            {
                size_t numOwners = other.owners.size();
                this->owners.resize(numOwners);
                for(size_t i = 0; i < numOwners; i++)
                {
                    this->owners[i] = other.owners[i];
                }        
                this->value = other.value;
                return (*this);
            };

            inline typename T::coord_type operator[](size_t idx) const{return this->value[idx];};

            inline typename T::coord_type &operator[](size_t idx){return this->value[idx];};

            inline OctTreeRankedValue operator+(const OctTreeRankedValue &other) const{return OctTreeRankedValue(this->value + other.value);};

            inline OctTreeRankedValue operator-(const OctTreeRankedValue &other) const{return OctTreeRankedValue(this->value - other.value);};

            inline OctTreeRankedValue operator*(double constant) const{return OctTreeRankedValue(this->value * constant);};

            inline OctTreeRankedValue operator/(double constant) const{return this->operator*(1/constant);};

            inline bool operator==(const OctTreeRankedValue &other) const{return this->value == other.value;};

            inline bool operator!=(const OctTreeRankedValue &other) const{return this->value != other.value;};

            friend std::ostream &operator<<(std::ostream &stream, const OctTreeRankedValue &wrapper)
            {
                size_t numOwners = wrapper.owners.size();
                if(numOwners == 0)
                {
                    return stream << wrapper.value << " [NO OWNER]";
                }
                if(numOwners == 1)
                {
                    return stream << wrapper.value << " [owner: " << wrapper.owners[0] << "]";
                }
                stream << wrapper.value << " [owners: ";
                for(size_t i = 0; i < numOwners - 1; i++)
                {
                    stream << wrapper.owners[i] << " ";
                }
                return stream << wrapper.owners[numOwners - 1] << "]";
            }
        };
    }
}

#endif // RANKED_VALUE_HPP