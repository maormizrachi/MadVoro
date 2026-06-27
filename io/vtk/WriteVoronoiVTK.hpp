#ifndef OUTPUT_VORONOI_WRITE_VTK_HPP
#define OUTPUT_VORONOI_WRITE_VTK_HPP

#ifdef MADVORO_WITH_VTK

#include <vector>
#include <string>
#include <filesystem>
#include "../../Voronoi3D.hpp"
#include "write_vtu_3d.hpp"

namespace fs = std::filesystem;

namespace MadVoro
{
  namespace IO
  {
    template <typename PointT>
    void WriteVoronoiVTK(const Voronoi3D<PointT> &tri, const std::string &filename, const std::vector<std::vector<double>> &data = std::vector<std::vector<double>>(), const std::vector<std::string>& names = std::vector<std::string>())
    {
        std::vector<std::vector<double>> vtu_cell_variables;
        std::vector<std::string> vtu_cell_variable_names;
        std::vector<std::string> vtu_cell_vectors_names;
        std::vector<std::vector<PointT>> vtu_cell_vectors;

        for(size_t i = 0; i < data.size(); ++i)
        {
            vtu_cell_variables.push_back(data[i]);
            vtu_cell_variable_names.push_back(names[i]);
        }

        size_t Npoints = tri.GetPointNo();
        vtu_cell_vectors_names.push_back("Coordinates");
        std::vector<PointT> vel(Npoints);
        for(size_t i = 0; i < Npoints; ++i)
            vel[i] = tri.GetMeshPoint(i);
        vtu_cell_vectors.push_back(vel);

        std::vector<double> temp(Npoints);
        for(size_t i = 0; i < Npoints; ++i)
            temp[i] = i;
        vtu_cell_variables.push_back(temp);
        vtu_cell_variable_names.push_back("Point Index");

        std::filesystem::path vtu_name(filename);
        vtu_name.replace_extension("vtu");
        MadVoro::IO::write_vtu3d::write_vtu_3d(vtu_name, vtu_cell_variable_names, vtu_cell_variables, vtu_cell_vectors_names, vtu_cell_vectors, tri);
    }
  }
}

#endif // MADVORO_WITH_VTK

#endif // OUTPUT_VORONOI_WRITE_VTK_HPP
