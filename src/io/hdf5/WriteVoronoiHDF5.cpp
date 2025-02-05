#include "WriteVoronoiHDF5.hpp"

#ifdef MADVORO_WITH_HDF5

namespace MadVoro
{
    namespace IO
    {
        void WriteVoronoiHDF5_Helper(H5File &file, Group &writegroup, const std::string &filename, const Voronoi3D &tri, const std::vector<std::vector<double>> &data, const std::vector<std::string> &names)
        {
            for(size_t i = 0; i < data.size(); ++i)
            {
                write_std_vector_to_hdf5(writegroup, data[i], names[i]);
            }

            vector<double> x, y, z, vx, vy, vz;
            vector<size_t> Nfaces;
            vector<size_t> Nvert;
            vector<size_t> FacesInCell;
            vector<size_t> VerticesInFace;
            size_t Npoints = tri.GetPointNo();

            for(size_t i = 0; i < Npoints; ++i)
            {
                const Vector3D &point = tri.GetMeshPoint(i);
                x.push_back(point.x);
                y.push_back(point.y);
                z.push_back(point.z);
            }

            write_std_vector_to_hdf5(writegroup, x, "mesh_point_x");
            write_std_vector_to_hdf5(writegroup, y, "mesh_point_y");
            write_std_vector_to_hdf5(writegroup, z, "mesh_point_z");
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
            write_std_vector_to_hdf5(writegroup, x, "all_mesh_point_x");
            write_std_vector_to_hdf5(writegroup, y, "all_mesh_point_y");
            write_std_vector_to_hdf5(writegroup, z, "all_mesh_point_z");

            for(size_t i = 0; i < Npoints; ++i)
            {
                const face_vec &face = tri.GetCellFaces(i);
                Nfaces.push_back(face.size());
                for(size_t j = 0; j < Nfaces.back(); ++j)
                {
                    FacesInCell.push_back(face[j]);
                }
            }
            H5::IntType datatype(PredType::NATIVE_ULLONG);
            datatype.setOrder(H5T_ORDER_LE);
            write_std_vector_to_hdf5(writegroup, Nfaces, "Number_of_faces_in_cell", datatype);
            write_std_vector_to_hdf5(writegroup, FacesInCell, "Faces_in_cell", datatype);
            Npoints = tri.GetFacePoints().size();
            for(size_t i = 0; i < Npoints; ++i)
            {
                vx.push_back(tri.GetFacePoints()[i].x);
                vy.push_back(tri.GetFacePoints()[i].y);
                vz.push_back(tri.GetFacePoints()[i].z);
            }
            write_std_vector_to_hdf5(writegroup, vx, "vertice_x");
            write_std_vector_to_hdf5(writegroup, vy, "vertice_y");
            write_std_vector_to_hdf5(writegroup, vz, "vertice_z");
            Npoints = tri.GetTotalFacesNumber();
            for(size_t i = 0; i < Npoints; ++i)
            {
                Nvert.push_back(tri.GetPointsInFace(i).size());
                for(size_t j = 0; j < Nvert.back(); ++j)
                {
                    VerticesInFace.push_back(tri.GetPointsInFace(i)[j]);
                }
            }
            write_std_vector_to_hdf5(writegroup, Nvert, "Number_of_vertices_in_face", datatype);
            write_std_vector_to_hdf5(writegroup, VerticesInFace, "Vertices_in_face", datatype);
        }

        #ifdef MADVORO_WITH_MPI
            void WriteVoronoiHDF5_Parallel(MadVoro::Voronoi3D const &tri, std::string const &filename, const std::vector<std::vector<double>> &data, const std::vector<std::string> &names)
            {
                int rank = 0;
                int ws = 0;
                H5File file;

                fs::path path = fs::absolute(filename).parent_path();
                std::string myFilePath;

                fs::path ranks_files_path = path / fs::path(filename).filename().replace_extension();
                if(not fs::exists(ranks_files_path))
                {
                    fs::create_directory(ranks_files_path);
                }
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                MPI_Comm_size(MPI_COMM_WORLD, &ws);
                myFilePath = (ranks_files_path / std::to_string(rank)).string() + ".h5";

                // truncate my file and open it
                H5File file2(H5std_string(myFilePath), H5F_ACC_TRUNC);
                file2.close();
                file.openFile(H5std_string(myFilePath), H5F_ACC_RDWR);

                Group writegroup = file.openGroup("/");

                if(rank == 0)
                {
                    std::vector<double> box(6);
                    box[0] = tri.GetBoxCoordinates().first.x;
                    box[1] = tri.GetBoxCoordinates().first.y;
                    box[2] = tri.GetBoxCoordinates().first.z;
                    box[3] = tri.GetBoxCoordinates().second.x;
                    box[4] = tri.GetBoxCoordinates().second.y;
                    box[5] = tri.GetBoxCoordinates().second.z;
                    write_std_vector_to_hdf5(writegroup, box, "Box");
                }

                WriteVoronoiHDF5_Helper(file, writegroup, filename, tri, data, names);

                writegroup.close();
                file.close();

                MPI_Barrier(MPI_COMM_WORLD);
                // only rank 0 makes the shared file
                if(rank == 0)
                {
                    file2 = H5File(H5std_string(filename), H5F_ACC_TRUNC);
                    file2.close();
                    hid_t shared_file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

                    for(int _rank = 0; _rank < ws; _rank++)
                    {
                        // merge `_rank`'s file
                        std::string rankFile((ranks_files_path / std::to_string(_rank)).string() + ".h5");
                        std::string rankGroupName("/rank" + std::to_string(_rank));
                        //Group rankGroup = sharedFile.createGroup(rankGroupName);
                        H5Lcreate_external(rankFile.c_str(),
                                        "/",
                                        shared_file_id,
                                        rankGroupName.c_str(),
                                        H5P_DEFAULT,
                                        H5P_DEFAULT);
                        //rankGroup.close();
                    }
                    H5Fclose(shared_file_id);
                }
            }
        #endif // MADVORO_WITH_MPI

        void WriteVoronoiHDF5(MadVoro::Voronoi3D const &tri, std::string const &filename, const std::vector<std::vector<double>> &data, const std::vector<std::string> &names)
        {
            #ifdef MADVORO_WITH_MPI
                int rank = 0;
                int ws = 0; // MPI_COMM_WORLD size
                MPI_Comm_rank(MPI_COMM_WORLD, &rank);
                MPI_Comm_size(MPI_COMM_WORLD, &ws);
            #endif

            H5File file;

            #ifdef MADVORO_WITH_MPI
                if(rank == 0)
                {
            #endif // MADVORO_WITH_MPI
                    H5File file2(H5std_string(filename), H5F_ACC_TRUNC);
                    file2.close();
                    file.openFile(H5std_string(filename), H5F_ACC_RDWR);
            #ifdef MADVORO_WITH_MPI
                }
            #endif // MADVORO_WITH_MPI

            Group writegroup;
            #ifdef MADVORO_WITH_MPI
                MPI_Barrier(MPI_COMM_WORLD);
                int dummy = 0;
                if(rank > 0)
                {
                    MPI_Recv(&dummy, 1, MPI_INT, rank - 1, HDF5_WRITE_BLOCK_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    file.openFile(H5std_string(filename), H5F_ACC_RDWR);
                }
                file.createGroup("/rank" + std::to_string(rank));
                writegroup = file.openGroup("/rank" + std::to_string(rank));
            #else
                writegroup = file.openGroup("/");
            #endif

            #ifdef MADVORO_WITH_MPI
                if(rank == 0)
                {
            #endif // MADVORO_WITH_MPI
                std::vector<double> box(6);
                box[0] = tri.GetBoxCoordinates().first.x;
                box[1] = tri.GetBoxCoordinates().first.y;
                box[2] = tri.GetBoxCoordinates().first.z;
                box[3] = tri.GetBoxCoordinates().second.x;
                box[4] = tri.GetBoxCoordinates().second.y;
                box[5] = tri.GetBoxCoordinates().second.z;
                write_std_vector_to_hdf5(writegroup, box, "Box");
            #ifdef MADVORO_WITH_MPI
                }
            #endif // MADVORO_WITH_MPI

            WriteVoronoiHDF5_Helper(file, writegroup, filename, tri, data, names);

            #ifdef MADVORO_WITH_MPI
                writegroup.close();
                file.close();
                if(rank < (ws - 1))
                {
                    int dummy = 0;
                    MPI_Send(&dummy, 1, MPI_INT, rank + 1, HDF5_WRITE_BLOCK_TAG, MPI_COMM_WORLD);
                }
                MPI_Barrier(MPI_COMM_WORLD);
            #else
                file.close();
            #endif
        }

        void WriteVoronoiHDF5_Serial(MadVoro::Voronoi3D const &tri, std::string const &filename, const std::vector<std::vector<double>> &data, const std::vector<std::string> &names)
        {
            H5File file(H5std_string(filename), H5F_ACC_TRUNC);
            std::vector<double> box(6);
            box[0] = tri.GetBoxCoordinates().first.x;
            box[1] = tri.GetBoxCoordinates().first.y;
            box[2] = tri.GetBoxCoordinates().first.z;
            box[3] = tri.GetBoxCoordinates().second.x;
            box[4] = tri.GetBoxCoordinates().second.y;
            box[5] = tri.GetBoxCoordinates().second.z;

            Group writegroup = file.openGroup("/");
            write_std_vector_to_hdf5(writegroup, box, "Box");
            WriteVoronoiHDF5_Helper(file, writegroup, filename, tri, data, names);
            writegroup.close();
            file.close();
        }
    }
}

#endif // MADVORO_WITH_HDF5