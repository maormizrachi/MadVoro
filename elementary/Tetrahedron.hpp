#ifndef MADVORO_TETRAHEDRON_HPP
#define MADVORO_TETRAHEDRON_HPP

#include <cstddef>

class Tetrahedron
{
public:
    std::size_t points[4];
    std::size_t neighbors[4];
    bool checkBig;
    bool newTetra;

    Tetrahedron() : points(), neighbors(), checkBig(true), newTetra(true) {}

    Tetrahedron(const Tetrahedron &other) = default;
    Tetrahedron &operator=(const Tetrahedron &other) = default;
    ~Tetrahedron() = default;
};

#endif // MADVORO_TETRAHEDRON_HPP
