#ifndef OUTPUT_VORONOI_WRITE_VTK_HPP
#define OUTPUT_VORONOI_WRITE_VTK_HPP

#ifdef MADVORO_WITH_VTK

#include <vector>
#include <string>
#include <filesystem>
#include "voronoi/Voronoi3D.hpp"
#include "write_vtu_3d.hpp"

namespace fs = std::filesystem;

namespace MadVoro
{
  namespace IO
  {
    /*! \brief Write voronoi data to a file, of all the ranks, but not parallely (excusively)
      \param tri Voronoit tessellation
      \param filename Name of output file
      \param write_vtu whether to write to vtu as well
    */
    void WriteVoronoiVTK(const MadVoro::Voronoi3D &tri, const std::string &filename, const std::vector<std::vector<double>> &data = std::vector<std::vector<double>>(), const std::vector<std::string>& names = std::vector<std::string>());
  }
}

#endif // MADVORO_WITH_VTK

#endif // OUTPUT_VORONOI_WRITE_VTK_HPP