/*! \file simple_io.hpp
  \brief A collection of simple input / output methods
  \author Almog Yalinewich
*/

#ifndef SIMPLE_IO_HPP
#define SIMPLE_IO_HPP 1

#include <iostream>
#include <fstream>

namespace MadVoro
{
  namespace IO
  {
    void binary_write_single_int(int n, std::ofstream &fh);

    /*! \brief Writes a double to a binary file
      \param d A double
      \param fh File handle
    */
    void binary_write_single_double(double d, std::ofstream &fh);

    /*! \brief Writes a single size_t to a binary file
      \param n size_t
      \param fh File handle
    */
    void binary_write_single_size_t(std::size_t n, std::ofstream &fh);
  }
}

#endif // SIMPLE_IO_HPP
