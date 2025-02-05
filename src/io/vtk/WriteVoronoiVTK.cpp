#include "WriteVoronoiVTK.hpp"

#ifdef MADVORO_WITH_VTK

using namespace MadVoro;

struct VTU_Output
{
    std::vector<std::vector<double>> vtu_cell_variables;
    std::vector<std::string> vtu_cell_variable_names;
    std::vector<std::string> vtu_cell_vectors_names;
    std::vector<std::vector<Vector3D>> vtu_cell_vectors;
};

void writeVTU(const std::string &filename, const MadVoro::Voronoi3D &tri, const VTU_Output &data)
{
    std::filesystem::path vtu_name(filename);
    vtu_name.replace_extension("vtu");
    MadVoro::IO::write_vtu3d::write_vtu_3d(vtu_name, data.vtu_cell_variable_names, data.vtu_cell_variables, data.vtu_cell_vectors_names, data.vtu_cell_vectors, tri);
}

VTU_Output WriteVoronoiVTKHelper(const std::string &filename, const MadVoro::Voronoi3D &tri, const std::vector<std::vector<double>> &data, const std::vector<std::string> &names)
{
    VTU_Output vtu;
    std::vector<std::vector<double>> &vtu_cell_variables = vtu.vtu_cell_variables;
    std::vector<std::string> &vtu_cell_variable_names = vtu.vtu_cell_variable_names;
    std::vector<std::string> &vtu_cell_vectors_names = vtu.vtu_cell_vectors_names;
    std::vector<std::vector<Vector3D>> &vtu_cell_vectors = vtu.vtu_cell_vectors;

    for(size_t i = 0; i < data.size(); ++i)
    {
        vtu_cell_variables.push_back(data[i]);
        vtu_cell_variable_names.push_back(names[i]);
    }

    std::vector<double> x, y, z, vx, vy, vz;
    std::vector<size_t> Nfaces;
    std::vector<size_t> Nvert;
    std::vector<size_t> FacesInCell;
    std::vector<size_t> VerticesInFace;
    size_t Npoints = tri.GetPointNo();

    vtu_cell_vectors_names.push_back("Coordinates");
    std::vector<Vector3D> vel(Npoints);
    for(size_t i = 0; i < Npoints; ++i)
        vel[i] = tri.GetMeshPoint(i);
    vtu_cell_vectors.push_back(vel);

    for(size_t i = 0; i < Npoints; ++i)
    {
        const Vector3D &point = tri.GetMeshPoint(i);
        x.push_back(point.x);
        y.push_back(point.y);
        z.push_back(point.z);
    }

    x.clear();
    y.clear();
    z.clear();

    for(size_t i = 0; i < tri.GetTotalPointNumber(); ++i)
    {
        const Vector3D &point = tri.GetMeshPoint(i);
        x.push_back(point.x);
        y.push_back(point.y);
        z.push_back(point.z);
    }

    std::vector<double> temp(Npoints);

    for(size_t i = 0; i < Npoints; ++i)
    {
        temp[i] = i;
    }
    vtu_cell_variables.push_back(temp);
    vtu_cell_variable_names.push_back("Point Index");

    for(size_t i = 0; i < Npoints; ++i)
    {
        const face_vec &face = tri.GetCellFaces(i);
        Nfaces.push_back(face.size());
        for(size_t j = 0; j < Nfaces.back(); ++j)
        {
            FacesInCell.push_back(face[j]);
        }
    }

    Npoints = tri.GetFacePoints().size();
    for(size_t i = 0; i < Npoints; ++i)
    {
        vx.push_back(tri.GetFacePoints()[i].x);
        vy.push_back(tri.GetFacePoints()[i].y);
        vz.push_back(tri.GetFacePoints()[i].z);
    }
    Npoints = tri.GetTotalFacesNumber();
    for(size_t i = 0; i < Npoints; ++i)
    {
        Nvert.push_back(tri.GetPointsInFace(i).size());
        for(size_t j = 0; j < Nvert.back(); ++j)
        {
            VerticesInFace.push_back(tri.GetPointsInFace(i)[j]);
        }
    }

    return vtu;
}

void MadVoro::IO::WriteVoronoiVTK(MadVoro::Voronoi3D const &tri, std::string const &filename, const std::vector<std::vector<double>> &data, const std::vector<std::string> &names)
{
    VTU_Output vtu = WriteVoronoiVTKHelper(filename, tri, data, names);
    writeVTU(filename, tri, vtu);
}

#endif // MADVORO_WITH_VTK