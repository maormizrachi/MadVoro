cmake_minimum_required(VERSION 3.14)
project(MadVoro
  VERSION 1.0.0
  LANGUAGES CXX
  DESCRIPTION "Massively distributed construction of 3D Voronoi diagrams"
)

# Only set these as defaults when MadVoro is the top-level project.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
      Debug Release RelWithDebInfo MinSizeRel)
  endif()

  set(CMAKE_CXX_STANDARD 17 CACHE STRING "")
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)

  # --- Compiler flags (aligned with RICH/config/compiler_flags.cmake) ---

  set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG -g1" CACHE STRING "" FORCE)
  set(CMAKE_CXX_FLAGS_DEBUG   "-O0 -g3 -gdwarf-3 -D_GLIBCXX_ASSERTIONS" CACHE STRING "" FORCE)

  string(APPEND CMAKE_CXX_FLAGS " -fno-omit-frame-pointer")
  string(APPEND CMAKE_CXX_FLAGS " -fstrict-aliasing")
  string(APPEND CMAKE_CXX_FLAGS " -fno-common")
  string(APPEND CMAKE_CXX_FLAGS " -rdynamic")
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
# into deps/.  The directory is auto-detected; override with -DMADVORO_DEPS_DIR=...
if(NOT DEFINED MADVORO_DEPS_DIR)
  set(MADVORO_DEPS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/deps")
endif()
get_filename_component(_deps "${MADVORO_DEPS_DIR}" ABSOLUTE)
if(IS_DIRECTORY "${_deps}")
  message(STATUS "MadVoro deps: ${_deps}")
  include_directories("${_deps}")
else()
  message(STATUS "MadVoro deps directory not found: ${_deps} (standalone deps unavailable)")
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
set(MADVORO_SOURCES "")

set(MADVORO_IO_SOURCES "")
if(MADVORO_WITH_HDF5)
  list(APPEND MADVORO_IO_SOURCES
    io/hdf5/hdf5_utils.cpp
  )
endif()

list(APPEND MADVORO_SOURCES
  ${MADVORO_IO_SOURCES}
  exception/MadVoroException.cpp
  exception/InvalidArgumentException.cpp
  exception/SizeException.cpp
)

if(MADVORO_WITH_MPI AND EXISTS "${_deps}/mpi_utils/AmountManager.cpp")
  list(APPEND MADVORO_SOURCES
    "${_deps}/mpi_utils/AmountManager.cpp"
  )
endif()

add_library(madvoro ${MADVORO_SOURCES})
add_library(MadVoro::madvoro ALIAS madvoro)

target_include_directories(madvoro PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# range_finders sub-submodule
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/range/finders")
  target_include_directories(madvoro PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/range/finders>
  )
endif()
target_link_libraries(madvoro PUBLIC Boost::boost OpenMP::OpenMP_CXX)

if(MADVORO_WITH_MPI)
  target_compile_definitions(madvoro PUBLIC MADVORO_WITH_MPI)
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
