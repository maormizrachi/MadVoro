#ifndef OUTPUT_POINTS_WRITE_HDF5_HPP
#define OUTPUT_POINTS_WRITE_HDF5_HPP

#ifdef MADVORO_WITH_HDF5

#include <string>
#include <filesystem>
#include "voronoi/Voronoi3D.hpp"
#include "hdf5_utils.hpp"

#ifdef MADVORO_WITH_MPI
  #include <mpi.h>
  #define HDF5_WRITE_BLOCK_TAG 604
#endif // MADVORO_WITH_MPI

namespace fs = std::filesystem;

using H5File = H5::H5File;

using namespace MadVoro;

namespace MadVoro
{
    namespace IO
    {
        void WritePointsHDF5(const std::vector<Vector3D> &points, const std::string &filename, const std::vector<std::vector<double>> &data = std::vector<std::vector<double>>(), const std::vector<std::string>& names = std::vector<std::string>());
    
        void WritePointsHDF5_Parallel(const std::vector<Vector3D> &points, const std::string &filename, const std::vector<std::vector<double>> &data = std::vector<std::vector<double>>(), const std::vector<std::string>& names = std::vector<std::string>());
    }
}

#endif // MADVORO_WITH_HDF5

#endif // OUTPUT_POINTS_WRITE_HDF5_HPP