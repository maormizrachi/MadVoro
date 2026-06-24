#ifndef OUTPUT_POINTS_WRITE_HDF5_HPP
#define OUTPUT_POINTS_WRITE_HDF5_HPP

#ifdef MADVORO_WITH_HDF5

#include <string>
#include <filesystem>
#include "voronoi/Voronoi3DFull.hpp"
#include "hdf5_utils.hpp"

#ifdef MADVORO_WITH_MPI
  #include <mpi.h>
  #define HDF5_WRITE_BLOCK_TAG 604
#endif // MADVORO_WITH_MPI

namespace fs = std::filesystem;

using H5File = H5::H5File;

namespace MadVoro
{
  namespace IO
  {
    #if MADVORO_WITH_MPI  
      void WriteVoronoiHDF5_Parallel(const MadVoro::Voronoi3DFull &tri, const std::string &filename, const std::vector<std::vector<double>> &data = std::vector<std::vector<double>>(), const std::vector<std::string>& names = std::vector<std::string>(), bool write_vtu = true);
    #endif // MADVORO_WITH_MPI

    void WriteVoronoiHDF5(const MadVoro::Voronoi3DFull &tri, const std::string &filename, const std::vector<std::vector<double>> &data = std::vector<std::vector<double>>(), const std::vector<std::string>& names = std::vector<std::string>(), bool write_vtu = true);

    void WriteVoronoiHDF5_Serial(const MadVoro::Voronoi3DFull &tri, const std::string &filename, const std::vector<std::vector<double>> &data = std::vector<std::vector<double>>(), const std::vector<std::string>& names = std::vector<std::string>(), bool write_vtu = true);
  }
}

#endif // MADVORO_WITH_HDF5

#endif // OUTPUT_POINTS_WRITE_HDF5_HPP
