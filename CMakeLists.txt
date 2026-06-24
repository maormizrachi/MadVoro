cmake_minimum_required(VERSION 3.14)
project(MadVoro
  VERSION 1.0.0
  LANGUAGES CXX
  DESCRIPTION "Massively distributed construction of 3D Voronoi diagrams"
)

# Only set these as defaults when MadVoro is the top-level project.
# When embedded via add_subdirectory(), the parent controls these.
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

# ---------- Required dependencies ----------
find_package(Boost 1.74.0 REQUIRED)

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

# ---------- Library ----------
add_subdirectory(src)

# ---------- Examples ----------
if(MADVORO_BUILD_EXAMPLES)
  add_subdirectory(examples)
endif()

# ---------- Install & export (skip when embedded via add_subdirectory) ----------
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

  install(DIRECTORY include/madvoro
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
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
