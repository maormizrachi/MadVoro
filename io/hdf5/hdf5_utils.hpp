/*! \file hdf5_utils.hpp
  \brief Higher level hdf5 utilities
  \author Almog Yalinewich
 */
#ifndef HDF5_UTILS_HPP
#define HDF5_UTILS_HPP 1

#ifdef MADVORO_WITH_HDF5

#include <string>
#include <vector>
#include <algorithm>
#include <H5Cpp.h>

using std::string;
using std::vector;
using std::pair;
using H5::Group;
using H5::PredType;
using H5::DataSpace;
using H5::DSetCreatPropList;
using H5::DataSet;
using H5::DataType;

namespace MadVoro
{
  namespace IO
  {
    /*! \brief Master function for writing vectors to hdf5 files
      \param file Either an actual file or a group within a file
      \param data Data to be written
      \param caption Name of dataset
      \param dt Data type
    */
    template<class T> void write_std_vector_to_hdf5(const Group& file, const vector<T>& data, const string& caption, const DataType& dt)
    {
      hsize_t dimsf[1];
      dimsf[0] = static_cast<hsize_t>(data.size());
      DataSpace dataspace(1, dimsf);

      DSetCreatPropList plist;
      if(dimsf[0]>100000)
        dimsf[0] = 100000;
      if (dimsf[0] == 0)
        dimsf[0] = 1;
      plist.setChunk(1,dimsf);
      plist.setDeflate(6);

      DataSet dataset = file.createDataSet
        (H5std_string(caption),
        dt,
        dataspace,
        plist);
      if(data.empty())
        dataset.write(nullptr, dt);
      else
      dataset.write(&data[0],dt);
    }

    /*! \brief Writes floating point data to hdf5
      \param file Either an actual file or a group within a file
      \param data Data to be written
      \param caption Name of dataset
    */
    void write_std_vector_to_hdf5(const Group& file, const vector<double>& data, const string& caption);

    /*! \brief Writes integer data to hdf5
      \param file Either an actual file or a group within a file
      \param data Data to be written
      \param caption Name of dataset
    */
    void write_std_vector_to_hdf5(const Group& file, const vector<int>& data, const string& caption);

    /*! \brief Writes size_t data to hdf5
    \param file Either an actual file or a group within a file
    \param data Data to be written
    \param caption Name of dataset
    */
    void write_std_vector_to_hdf5(const Group& file, const vector<size_t>& data, const string& caption);

    //! \brief Facilitates writing hdf5 files
    class HDF5Shortcut
    {
    public:

      /*! \brief Class constructor
        \param fname Name of hdf5 file
      */
      explicit HDF5Shortcut(const string& fname);

      /*! \brief adds dataset
        \param field_name Name of dataset
        \param data_set Array of data
        \return Self reference
      */
      HDF5Shortcut& operator()(const string& field_name,
            const vector<double>& data_set);
      
      /*! \brief adds dataset
        \param field_name Name of dataset
        \param data_set Array of data
        \return Self reference
      */
      HDF5Shortcut& operator()(const string &field_name, const vector<int> &data_set);

      //! \brief Class destructor. This is the stage when the file is written
      ~HDF5Shortcut(void);

    private:
      const string fname_;
      vector<pair<string,vector<double> > > double_data_;
      vector<pair<string,vector<int> > > int_data_;
    };
  }
}

#endif // MADVORO_WITH_HDF5

#endif // HDF5_UTILS_HPP
