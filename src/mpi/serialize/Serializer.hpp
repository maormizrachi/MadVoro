#ifndef SERIALIZER_HPP
#define SERIALIZER_HPP

#ifdef MADVORO_WITH_MPI

#include <array>
#include <vector>
#include <memory>
#include <limits>
#include <cstring>
#include <cassert>

#include "Serializable.hpp"

#include "utils/compiler.h"

namespace MadVoro
{
    namespace MPI
    {
        template<typename Test, template<typename...> class Ref>
        struct is_specialization : std::false_type {};

        template<template<typename...> class Ref, typename... Args>
        struct is_specialization<Ref<Args...>, Ref>: std::true_type {};

        template<typename T>
        using is_vector = is_specialization<T, std::vector>;

        class Serializer
        {
        public:
            Serializer(void) = default;

            void reset(){this->internal.clear();};

            char *resize(size_t size);

            template<typename T>
            size_t insert(const T &data);

            template<typename T, template<typename...> class VectorContainer, typename... Ts>
            size_t insert(const VectorContainer<T, Ts...> &data, size_t startIndex, size_t writeSize);

            template<typename T, template<typename...> class VectorContainer, typename... Ts>
            inline size_t insert_all(const VectorContainer<T, Ts...> &data, size_t offset = 0){return this->insert(data, offset, std::numeric_limits<size_t>::max());};

            template<typename T>
            size_t insert_array(const T *data, size_t arraySize);

            template<typename T, size_t N>
            inline size_t insert_array(const std::array<T, N> &data){return this->insert_array(data.data(), N);};
            
            template<typename T, typename Index_T = size_t>
            size_t insert_all_indexed(const std::vector<T> &data, const std::vector<Index_T> &indices, size_t extent = 1);
            
            template<typename T>
            size_t insert_elements(const std::vector<T> &data, size_t startIndex, size_t numElements);

            template<typename T>
            size_t extract(T &data, size_t idx) const;

            template<typename T, template<typename...> class VectorContainer, typename... Ts>
            size_t extract(VectorContainer<T, Ts...> &values, size_t startIndex, size_t endIndex) const;

            template<typename T, template<typename...> class VectorContainer, typename... Ts>
            inline size_t extract_all(VectorContainer<T, Ts...> &values, size_t offset = 0) const{return this->extract(values, offset, this->internal.size());};

            template<typename T>
            size_t extract_array(T *values, size_t arraySize, size_t offset = 0) const;

            template<typename T, size_t N>
            inline size_t extract_array(std::array<T, N> &values, size_t offset = 0) const{return this->extract_array(values.data(), N, offset);};

            inline const char *getData() const{return this->internal.data();};

            inline char *getData(){return this->internal.data();};

            inline size_t size() const{return this->internal.size();};

        private:
            std::vector<char> internal; // internal buffer
        };

        inline char *Serializer::resize(size_t size)
        {
            size_t oldSize = this->internal.size();
            this->internal.resize(oldSize + size);
            return this->internal.data() + oldSize;
        }

        template<typename T>
        force_inline size_t Serializer::insert(const T &data)
        {
            if constexpr(is_serializable<T>::value)
            {
                // serializable
                return data.dump(this);
            }
            else
            {
                if constexpr(is_vector<T>::value)
                {
                    // vector
                    size_t bytes = 0;
                    using S = typename T::value_type;
                    const std::vector<S> &vector = reinterpret_cast<const std::vector<S>&>(data);
                    size_t putSizeIdx = this->internal.size();
                    bytes += this->insert<size_t>(static_cast<size_t>(0));
                    size_t vectorBytes = 0;
                    vectorBytes = this->insert_all(vector);
                    bytes += vectorBytes;
                    *(reinterpret_cast<size_t*>(this->internal.data() + putSizeIdx)) = vectorBytes;
                    return bytes;
                }
                else
                {
                    // not serializable
                    constexpr size_t size = sizeof(T);
                    const char *ptr = reinterpret_cast<const char*>(&data);
                    this->internal.insert(this->internal.end(), ptr, ptr + size);
                    return size;
                }
            }
        }

        template<typename T, template<typename...> class VectorContainer, typename... Ts>
        force_inline size_t Serializer::insert(const VectorContainer<T, Ts...> &data, size_t startIndex, size_t writeSize)
        {
            size_t index = startIndex;
            size_t bytes = 0;
            while((bytes < writeSize) and (index < data.size()))
            {
                bytes += this->insert(data[index]);
                index++;
            }
            return bytes;
        }

        template<typename T>
        force_inline size_t Serializer::insert_elements(const std::vector<T> &data, size_t startIndex, size_t numElements)
        {
            size_t index = startIndex;
            size_t inserted = 0;
            size_t bytes = 0;
            while((inserted < numElements) and (index < data.size()))
            {
                bytes += this->insert(data[index]);
                index++;
                inserted++;
            }
            return bytes;
        }

        template<typename T, typename Index_T>
        force_inline size_t Serializer::insert_all_indexed(const std::vector<T> &data, const std::vector<Index_T> &indices, size_t extent)
        {
            size_t bytes = 0;
            size_t N = indices.size();
            for(size_t i = 0; i < N; ++i)
            {
                for(size_t j = 0; j < extent; ++j)
                    bytes += this->insert(data[indices[i] * extent + j]);
            }
            return bytes;
        }

        template<typename T>
        force_inline size_t Serializer::insert_array(const T *data, size_t arraySize)
        {
            size_t bytes = 0;
            for(size_t i = 0; i < arraySize; ++i)
            {
                bytes += this->insert(data[i]);
            }
            return bytes;
        }

        template<typename T>
        force_inline size_t Serializer::extract(T &data, size_t idx) const
        {
            assert(idx < this->internal.size());
            if constexpr(is_serializable<T>::value)
            {
                // serializable
                return data.load(this, idx);
            }
            else
            {
                if constexpr(is_vector<T>::value)
                {
                    // vector
                    size_t bytes = 0;
                    using S = typename T::value_type;
                    std::vector<S> &vector = reinterpret_cast<std::vector<S>&>(data);
                    vector.clear();
                    size_t size;
                    bytes += this->extract(size, idx);
                    bytes += this->extract(vector, idx + bytes, size);
                    return bytes;
                }
                else
                {
                    // not serializable
                    constexpr size_t size = sizeof(T);
                    char *ptr = reinterpret_cast<char*>(&data);
                    std::memcpy(ptr, this->internal.data() + idx, size);
                    return size;
                }
            }
        }

        template<typename T, template<typename...> class VectorContainer, typename... Ts>
        force_inline size_t Serializer::extract(VectorContainer<T, Ts...> &values, size_t startOffset, size_t readSize) const
        {
            size_t bytesRead = 0;
            while((bytesRead < readSize) and ((startOffset + bytesRead) < this->internal.size()))
            {
                values.emplace_back();
                bytesRead += this->extract(values.back(), startOffset + bytesRead);
            }
            return bytesRead;
        }

        template<typename T>
        force_inline size_t Serializer::extract_array(T *values, size_t arraySize, size_t offset) const
        {
            size_t bytesRead = 0;
            size_t index = 0;
            while(((offset + bytesRead) < this->internal.size()) and (index < arraySize))
            {
                bytesRead += this->extract(values[index], offset + bytesRead);
                index++;
            }
            return bytesRead;
        }
    }
}

#endif // MADVORO_WITH_MPI

#endif // SERIALIZER_HPP