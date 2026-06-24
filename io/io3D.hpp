/*! \file io3D.hpp
\brief A collection of simple input / output methods for 3D vectors
\author Elad Steinberg
*/

#ifndef IO3D_HPP
#define IO3D_HPP 1

#include "simple_io.hpp"
#include "3D/elementary/Vector3D.hpp"

namespace MadVoro
{
    namespace IO
    {
        void write_vec3d(std::vector<Vector3D> const&vec, std::string const& fname);

        std::vector<Vector3D> read_vec3d(std::string fname);

        void write_vecst(std::vector<size_t> const&vec, std::string const& fname);

        void write_vecint(std::vector<int> const&vec, std::string const& fname);

        std::vector<size_t> read_vecst(std::string fname);

        std::vector<int> read_vecint(std::string fname);
    }
}

#endif // IO3D_HPP
