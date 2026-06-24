cmake_minimum_required(VERSION 3.14)
project(MadVoro
  VERSION 1.0.0
  LANGUAGES CXX
  DESCRIPTION "Massively distributed construction of 3D Voronoi diagrams"
)

# Only set these as defaults when MadVoro is the top-level project.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  set(CMAKE_CXX_STANDARD 17 CACHE STRING "")
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
endif()

# ---------- Options ----------
option(MADVORO_WITH_MPI   "Build with MPI support"   OFF)
option(MADVORO_WITH_HDF5  "Build with HDF5 support"  OFF)
option(MADVORO_WITH_VTK   "Build with VTK support"   OFF)
option(MADVORO_WITH_VCL   "Build with VCL (Vector Class Library) support" OFF)
option(MADVORO_BUILD_EXAMPLES "Build example programs" OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

# ---------- External dependencies directory ----------
# When building standalone, users run install_deps.sh to clone dependencies
# into deps/, then pass -DMADVORO_DEPS_DIR=deps to CMake.
if(DEFINED MADVORO_DEPS_DIR)
  get_filename_component(_deps "${MADVORO_DEPS_DIR}" ABSOLUTE)
  include_directories("${_deps}/mpi_utils")
  include_directories("${_deps}/spatial_ds")
  include_directories("${_deps}/MeshDecomposer3D")
endif()

# ---------- Required dependencies ----------
find_package(Boost 1.74.0 REQUIRED)
find_package(OpenMP REQUIRED)

# ---------- Optional dependencies ----------
if(MADVORO_WITH_MPI)
  find_package(MPI 2.0 REQUIRED COMPONENTS CXX)
endif()

if(MADVORO_WITH_HDF5)
  set(HDF5_PREFER_PARALLEL ${MADVORO_WITH_MPI})
  find_package(HDF5 1.8.0 REQUIRED COMPONENTS CXX HL)
endif()

if(MADVORO_WITH_VTK)
  find_package(VTK 9.3.0 REQUIRED
    COMPONENTS
      CommonCore CommonColor FiltersGeneral FiltersSources
      IOXML InteractionStyle
  )
  if(MADVORO_WITH_MPI)
    find_package(VTK REQUIRED COMPONENTS IOParallelXML ParallelMPI)
  endif()
endif()

if(MADVORO_WITH_VCL)
  find_package(VCL REQUIRED)
endif()

# ---------- Library sources ----------
set(MADVORO_SOURCES
  Voronoi3D.cpp
  delaunay/Delaunay3D.cpp
  geometry/Intersections.cpp
  utils/Predicates3D.cpp
  utils/Predicates3D_old.cpp
)

set(MADVORO_IO_SOURCES "")
if(MADVORO_WITH_HDF5)
  list(APPEND MADVORO_IO_SOURCES
    io/hdf5/hdf5_utils.cpp
    io/hdf5/WriteVoronoiHDF5.cpp
    io/hdf5/WritePointsHDF5.cpp
  )
endif()
if(MADVORO_WITH_VTK)
  list(APPEND MADVORO_IO_SOURCES
    io/vtk/write_vtu_3d.cpp
    io/vtk/WriteVoronoiVTK.cpp
  )
endif()

list(APPEND MADVORO_SOURCES
  io/simple_io.cpp
  io/io3D.cpp
  ${MADVORO_IO_SOURCES}
  exception/MadVoroException.cpp
  exception/InvalidArgumentException.cpp
  exception/SizeException.cpp
)

add_library(madvoro ${MADVORO_SOURCES})
target_include_directories(madvoro PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)
target_link_libraries(madvoro PUBLIC Boost::boost OpenMP::OpenMP_CXX)

if(MADVORO_WITH_MPI)
  target_compile_definitions(madvoro PUBLIC MADVORO_WITH_MPI RICH_MPI)
  target_link_libraries(madvoro PUBLIC MPI::MPI_CXX)
endif()

if(MADVORO_WITH_HDF5)
  target_compile_definitions(madvoro PUBLIC MADVORO_WITH_HDF5)
  target_include_directories(madvoro PUBLIC ${HDF5_INCLUDE_DIRS})
  target_link_libraries(madvoro PUBLIC ${HDF5_LIBRARIES})
endif()

if(MADVORO_WITH_VTK)
  target_compile_definitions(madvoro PUBLIC MADVORO_WITH_VTK)
  target_link_libraries(madvoro PUBLIC ${VTK_LIBRARIES})
endif()

if(MADVORO_WITH_VCL)
  target_compile_definitions(madvoro PUBLIC USE_VCL_VECTORIZATION)
  target_include_directories(madvoro PUBLIC ${VCL_INCLUDE_DIR})
endif()

# ---------- Examples ----------
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR AND MADVORO_BUILD_EXAMPLES)
  include(CTest)
endif()

if(MADVORO_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

# ---------- Install & export (skip when embedded) ----------
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  include(GNUInstallDirs)
  include(CMakePackageConfigHelpers)

  install(TARGETS madvoro
    EXPORT MadVoroTargets
    LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME  DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  )

  install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/madvoro
    FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h"
    PATTERN ".git" EXCLUDE
    PATTERN "examples" EXCLUDE
    PATTERN "cmake" EXCLUDE
  )

  install(EXPORT MadVoroTargets
    FILE MadVoroTargets.cmake
    NAMESPACE MadVoro::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MadVoro
  )

  configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/MadVoroConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/MadVoroConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MadVoro
  )

  write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/MadVoroConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
  )

  install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/MadVoroConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/MadVoroConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/MadVoro
  )
endif()
