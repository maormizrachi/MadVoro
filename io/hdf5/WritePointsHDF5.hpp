#ifndef OUTPUT_POINTS_WRITE_POINTS_HDF5_HPP
#define OUTPUT_POINTS_WRITE_POINTS_HDF5_HPP

#ifdef MADVORO_WITH_HDF5

#include <string>
#include <filesystem>
#include "hdf5_utils.hpp"

#ifdef MADVORO_WITH_MPI
  #include <mpi.h>
  #define HDF5_WRITE_POINTS_BLOCK_TAG 605
#endif // MADVORO_WITH_MPI

namespace fs = std::filesystem;

using H5File = H5::H5File;

namespace MadVoro
{
    namespace IO
    {
        template <typename PointT>
        void WritePointsHDF5(const std::vector<PointT> &points, const std::string &filename, const std::vector<std::vector<double>> &data = std::vector<std::vector<double>>(), const std::vector<std::string>& names = std::vector<std::string>())
        {
            int rank = 0;
            int ws = 0;
            H5File file;
        #ifdef MADVORO_WITH_MPI
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &ws);
            if(rank == 0)
            {
        #endif
                H5File file2(H5std_string(filename), H5F_ACC_TRUNC);
                file2.close();
                file.openFile(H5std_string(filename), H5F_ACC_RDWR);
        #ifdef MADVORO_WITH_MPI
            }
        #endif
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

            size_t const Ncells = points.size();
            vector<double> temp(Ncells);
            for(size_t i = 0; i < Ncells; ++i)
                temp[i] = points[i].x;
            write_std_vector_to_hdf5(writegroup, temp, "X");
            for(size_t i = 0; i < Ncells; ++i)
                temp[i] = points[i].y;
            write_std_vector_to_hdf5(writegroup, temp, "Y");
            for(size_t i = 0; i < Ncells; ++i)
                temp[i] = points[i].z;
            write_std_vector_to_hdf5(writegroup, temp, "Z");
            for(size_t i = 0; i < data.size(); ++i)
                write_std_vector_to_hdf5(writegroup, data[i], names[i]);
            for(size_t i = 0; i < Ncells; ++i)
                temp[i] = rank;
            write_std_vector_to_hdf5(writegroup, temp, "MPI_rank");

        #ifdef MADVORO_WITH_MPI
            if(rank < (ws - 1))
            {
                int dummy = 0;
                writegroup.close();
                file.close();
                MPI_Send(&dummy, 1, MPI_INT, rank + 1, HDF5_WRITE_BLOCK_TAG, MPI_COMM_WORLD);
            }
            else
            {
                writegroup.close();
                file.close();
            }
            MPI_Barrier(MPI_COMM_WORLD);
        #else
            writegroup.close();
            file.close();
        #endif
        }

#ifdef MADVORO_WITH_MPI
        template <typename PointT>
        void WritePointsHDF5_Parallel(const std::vector<PointT> &points, const std::string &filename, const std::vector<std::vector<double>> &data = std::vector<std::vector<double>>(), const std::vector<std::string>& names = std::vector<std::string>())
        {
            int rank = 0;
            int ws = 0;
            H5File file;

            fs::path path = fs::absolute(filename).parent_path();
            std::string myFilePath;
            fs::path ranks_files_path = path / fs::path(filename).filename().replace_extension();
            if(not fs::exists(ranks_files_path))
                fs::create_directory(ranks_files_path);
            MPI_Comm_rank(MPI_COMM_WORLD, &rank);
            MPI_Comm_size(MPI_COMM_WORLD, &ws);
            myFilePath = (ranks_files_path / std::to_string(rank)).string() + ".h5";

            H5File file2(H5std_string(myFilePath), H5F_ACC_TRUNC);
            file2.close();
            file.openFile(H5std_string(myFilePath), H5F_ACC_RDWR);
            Group writegroup = file.openGroup("/");

            const size_t Ncells = points.size();
            vector<double> temp(Ncells);
            for(size_t i = 0; i < Ncells; ++i)
                temp[i] = points[i].x;
            write_std_vector_to_hdf5(writegroup, temp, "X");
            for(size_t i = 0; i < Ncells; ++i)
                temp[i] = points[i].y;
            write_std_vector_to_hdf5(writegroup, temp, "Y");
            for(size_t i = 0; i < Ncells; ++i)
                temp[i] = points[i].z;
            write_std_vector_to_hdf5(writegroup, temp, "Z");
            for(size_t i = 0; i < data.size(); ++i)
                write_std_vector_to_hdf5(writegroup, data[i], names[i]);
            for(size_t i = 0; i < Ncells; ++i)
                temp[i] = rank;
            write_std_vector_to_hdf5(writegroup, temp, "MPI_rank");

            writegroup.close();
            file.close();

            MPI_Barrier(MPI_COMM_WORLD);
            if(rank == 0)
            {
                file2 = H5File(H5std_string(filename), H5F_ACC_TRUNC);
                file2.close();
                hid_t shared_file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
                for(int _rank = 0; _rank < ws; _rank++)
                {
                    std::string rankFile((ranks_files_path / std::to_string(_rank)).string() + ".h5");
                    std::string rankGroupName("/rank" + std::to_string(_rank));
                    H5Lcreate_external(rankFile.c_str(), "/", shared_file_id, rankGroupName.c_str(), H5P_DEFAULT, H5P_DEFAULT);
                }
                H5Fclose(shared_file_id);
            }
        }
#endif // MADVORO_WITH_MPI
    }
}

#endif // MADVORO_WITH_HDF5

#endif // OUTPUT_POINTS_WRITE_POINTS_HDF5_HPP
