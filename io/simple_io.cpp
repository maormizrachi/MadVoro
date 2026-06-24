#include "io/simple_io.hpp"
#include "exception/MadVoroException.hpp"
#include <cassert>

void MadVoro::IO::binary_write_single_int(int n, std::ofstream &fh)
{
  fh.write(reinterpret_cast<const char*>(&n),sizeof(int));
}

void MadVoro::IO::binary_write_single_double(double d, std::ofstream &fh)
{
  fh.write(reinterpret_cast<const char*>(&d),sizeof(double));
}

void MadVoro::IO::binary_write_single_size_t(size_t n, std::ofstream &fh)
{
  fh.write(reinterpret_cast<const char*>(&n),sizeof(size_t));
}
