#ifndef PRINT_VECTORS_HPP
#define PRINT_VECTORS_HPP

#include <ostream>
#include <vector>
#include <boost/container/small_vector.hpp>

#include "printTuples.hpp"

#define ELEMENTS_TO_PRINT_VECTOR 3

namespace MadVoro
{
    template<typename ContainerType>
    std::ostream &printVec(std::ostream &stream, const ContainerType &vector)
    {
        if(vector.empty())
        {
            return stream << "{}";
        }

        stream << "{";
        size_t firstElementsToShow = std::min<size_t>(vector.size(), ELEMENTS_TO_PRINT_VECTOR);
        for(size_t i = 0; i < firstElementsToShow - 1; i++)
        {
            stream << vector[i] << ", ";
        }
        stream << vector[firstElementsToShow - 1];
        if(firstElementsToShow < vector.size())
        {
            if(firstElementsToShow == vector.size() - 1)
            {
                stream << ", ";
            }
            else
            {
                stream << ", ... ,";
            }
            size_t lastElementsToShow = std::min<size_t>(vector.size() - firstElementsToShow, ELEMENTS_TO_PRINT_VECTOR);
            for(size_t i = vector.size() - lastElementsToShow; i < vector.size() - 1; i++)
            {
                stream << vector[i] << ", ";
            }
            stream << vector[vector.size() - 1];
        }
        return stream << "}";
    }

    template<typename T>
    std::ostream &operator<<(std::ostream &stream, const std::vector<T> &vector)
    {
        return printVec(stream, vector);
    }

    template<typename T, size_t N, typename... Args>
    std::ostream &operator<<(std::ostream &stream, const boost::container::small_vector<T, N, Args...> &vector)
    {
        return printVec(stream, vector);
    }
}

#endif // PRINT_VECTORS_HPP