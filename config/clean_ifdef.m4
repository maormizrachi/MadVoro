AC_CONFIG_COMMANDS([clean_ifdef_directives], [
    # Always copy the file first
    mkdir -p src/include/madvoro
    cp -f src/elementary/Face.hpp src/include/madvoro/Face.hpp
    cp -f src/elementary/Vector3D.hpp src/include/madvoro/Vector3D.hpp
    cp -f src/voronoi/Voronoi3D.hpp src/include/madvoro/Voronoi3D.hpp

    if [[ "$(uname)" == "Darwin" ]]; then
        SED_CMD="sed -i ''"  # macOS needs an additional parameter in case of using the -i flag.
    else  
        SED_CMD="sed -i"
    fi

    $SED_CMD 's|elementary/||g' src/include/madvoro/Voronoi3D.hpp

    if test "$mpi_enabled" = "yes"; then
        # If mpi_enabled is yes, remove the `#ifdef` and `#endif` but keep the content inside
		$SED_CMD '/#ifdef MADVORO_WITH_MPI/d' src/include/madvoro/Voronoi3D.hpp
		$SED_CMD '/#endif \/\/ MADVORO_WITH_MPI/d' src/include/madvoro/Voronoi3D.hpp
	else
        # If mpi_enabled is no, remove everything inside the `#ifdef` block
		$SED_CMD '/#ifdef MADVORO_WITH_MPI/,/#endif \/\/ MADVORO_WITH_MPI/d' src/include/madvoro/Voronoi3D.hpp
    fi
    if test "$hdf5_enabled" = "yes"; then
		$SED_CMD '/#ifdef MADVORO_WITH_HDF5/d' src/include/madvoro/Voronoi3D.hpp
		$SED_CMD '/#endif \/\/ MADVORO_WITH_HDF5/d' src/include/madvoro/Voronoi3D.hpp
	else
		$SED_CMD '/#ifdef MADVORO_WITH_HDF5/,/#endif \/\/ MADVORO_WITH_HDF5/d' src/include/madvoro/Voronoi3D.hpp
    fi
    if test "$vtk_enabled" = "yes"; then
		$SED_CMD '/#ifdef MADVORO_WITH_VTK/d' src/include/madvoro/Voronoi3D.hpp
		$SED_CMD '/#endif \/\/ MADVORO_WITH_VTK/d' src/include/madvoro/Voronoi3D.hpp
	else
		$SED_CMD '/#ifdef MADVORO_WITH_VTK/,/#endif \/\/ MADVORO_WITH_VTK/d' src/include/madvoro/Voronoi3D.hpp
    fi
], [mpi_enabled=$mpi_enabled hdf5_enabled=$hdf5_enabled vtk_enabled=$vtk_enabled])
