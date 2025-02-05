#ifndef PRINT_TUPLES_HPP
#define PRINT_TUPLES_HPP

#include <ostream>
#include <tuple>

namespace MadVoro
{
    template<typename... Args>
    std::ostream &operator<<(std::ostream &stream, const std::tuple<Args...> &tuple)
    {
        stream << "(";
        std::apply([&stream](auto &&... args){((stream << args << ", "), ...);}, tuple);
        return stream << ")";
    }

    template<typename T, typename U>
    std::ostream &operator<<(std::ostream &stream, const std::pair<T, U> &pair)
    {
        return stream << "(" << pair.first << ", " << pair.second << ")";
    }
}

#endif // PRINT_TUPLES_HPP