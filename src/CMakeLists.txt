add_library(madvoro
  delaunay/Delaunay3D.cpp
  elementary/Face3D.cpp
  elementary/Point3D.cpp
  elementary/Tetrahedron.cpp
  exception/InvalidArgumentException.cpp
  exception/MadVoroException.cpp
  exception/SizeException.cpp
  geometry/Intersections.cpp
  hilbert/HilbertOrder3D.cpp
  hilbert/HilbertOrder3D_Utils.cpp
  hilbert/rectangular/HilbertRectangularConvertor3D.cpp
  io/io3D.cpp
  io/simple_io.cpp
  io/hdf5/hdf5_utils.cpp
  io/hdf5/WritePointsHDF5.cpp
  io/hdf5/WriteVoronoiHDF5.cpp
  io/vtk/write_vtu_3d.cpp
  io/vtk/WriteVoronoiVTK.cpp
  loadBalancing/CurveLoadBalancer.cpp
  loadBalancing/HilbertLoadBalancer.cpp
  loadBalancing/LoadBalancer.cpp
  utils/Predicates3D.cpp
  utils/utils.cpp
  voronoi/Voronoi3DFull.cpp
  voronoi/pointsManager/HilbertPointsManager.cpp
  voronoi/pointsManager/PointsManager.cpp
)

add_library(MadVoro::madvoro ALIAS madvoro)

target_include_directories(madvoro
  PUBLIC
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:include>
)

target_compile_features(madvoro PUBLIC cxx_std_17)

target_link_libraries(madvoro PUBLIC Boost::boost)

# ---------- MPI ----------
if(MADVORO_WITH_MPI)
  target_compile_definitions(madvoro PUBLIC MADVORO_WITH_MPI)
  target_link_libraries(madvoro PUBLIC MPI::MPI_CXX)
endif()

# ---------- HDF5 ----------
if(MADVORO_WITH_HDF5)
  target_compile_definitions(madvoro PUBLIC MADVORO_WITH_HDF5)
  target_include_directories(madvoro PRIVATE ${HDF5_INCLUDE_DIRS})
  target_link_libraries(madvoro PRIVATE ${HDF5_CXX_LIBRARIES} ${HDF5_HL_LIBRARIES} ${HDF5_LIBRARIES})
endif()

# ---------- VTK ----------
if(MADVORO_WITH_VTK)
  target_compile_definitions(madvoro PUBLIC MADVORO_WITH_VTK)
  target_link_libraries(madvoro PRIVATE ${VTK_LIBRARIES})
  vtk_module_autoinit(TARGETS madvoro MODULES ${VTK_LIBRARIES})
endif()

# ---------- VCL ----------
if(MADVORO_WITH_VCL)
  target_compile_definitions(madvoro PUBLIC MADVORO_WITH_VCL)
  target_include_directories(madvoro PRIVATE ${VCL_INCLUDE_DIRS})
endif()

set_target_properties(madvoro PROPERTIES
  VERSION   ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
  POSITION_INDEPENDENT_CODE ON
)
