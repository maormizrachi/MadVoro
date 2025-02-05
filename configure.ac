AC_INIT([madvoro],[1.0])

AC_CONFIG_SRCDIR([src])

AM_INIT_AUTOMAKE([1.9 foreign -Wall, -Werror, subdir-objects])

AM_SILENT_RULES([yes])

AM_PROG_AR
AC_PROG_CXX
LT_INIT([shared], [4.1.6])  # Ensure shared library support

m4_include([config/check_mpi.m4])
m4_include([config/check_hdf5.m4])
m4_include([config/check_boost.m4])
m4_include([config/check_vcl.m4])
m4_include([config/check_vtk.m4])

AC_LANG_PUSH([C++])  # Switch to C++

# Handle MPI

MPI_HANDLE([2.0], [mpi_enabled=yes;
                    AC_DEFINE([MADVORO_WITH_MPI], [1], [Define if building with MPI support])
                    ], [mpi_enabled=no])

BOOST_HANDLE([1.74.0], [boost_found=yes], [boost_found=no])
if test "$boost_found" = "no"; then
    AC_MSG_ERROR("Boost was not found")
fi

HDF5_HANDLE([1.8.0], [hdf5_enabled=yes; 
                        AC_DEFINE([MADVORO_WITH_HDF5], [1], [Define if building with HDF5 support])
                     ], [hdf5_enabled=no])

VCL_HANDLE([20000], [vcl_enabled=yes;
                        AC_DEFINE([MADVORO_WITH_VCL], [1], [Define if building with VCL support])
                    ], [vcl_enabled=no])

VTK_HANDLE([9.3.0], [$mpi_enabled], [vtk_enabled=yes; 
                                        AC_DEFINE([MADVORO_WITH_VTK], [1], [Define if building with VTK support])
                                    ], [vtk_enabled=no])

AC_CONFIG_COMMANDS([clean_ifdef_directives], [
    # Always copy the file first
    mkdir -p src/include/madvoro
    cp -f src/elementary/Face.hpp src/include/madvoro/Face.hpp
    cp -f src/elementary/Vector3D.hpp src/include/madvoro/Vector3D.hpp
    cp -f src/voronoi/Voronoi3D.hpp src/include/madvoro/Voronoi3D.hpp
    sed -i 's|elementary/||g' src/include/madvoro/Voronoi3D.hpp

    if test "$mpi_enabled" = "yes"; then
        # If mpi_enabled is yes, remove the `#ifdef` and `#endif` but keep the content inside
		sed -i '/#ifdef MADVORO_WITH_MPI/d' src/include/madvoro/Voronoi3D.hpp
		sed -i '/#endif \/\/ MADVORO_WITH_MPI/d' src/include/madvoro/Voronoi3D.hpp
	else
        # If mpi_enabled is no, remove everything inside the `#ifdef` block
		sed -i '/#ifdef MADVORO_WITH_MPI/,/#endif \/\/ MADVORO_WITH_MPI/d' src/include/madvoro/Voronoi3D.hpp
    fi
    if test "$hdf5_enabled" = "yes"; then
		sed -i '/#ifdef MADVORO_WITH_HDF5/d' src/include/madvoro/Voronoi3D.hpp
		sed -i '/#endif \/\/ MADVORO_WITH_HDF5/d' src/include/madvoro/Voronoi3D.hpp
	else
		sed -i '/#ifdef MADVORO_WITH_HDF5/,/#endif \/\/ MADVORO_WITH_HDF5/d' src/include/madvoro/Voronoi3D.hpp
    fi
    if test "$vtk_enabled" = "yes"; then
		sed -i '/#ifdef MADVORO_WITH_VTK/d' src/include/madvoro/Voronoi3D.hpp
		sed -i '/#endif \/\/ MADVORO_WITH_VTK/d' src/include/madvoro/Voronoi3D.hpp
	else
		sed -i '/#ifdef MADVORO_WITH_VTK/,/#endif \/\/ MADVORO_WITH_VTK/d' src/include/madvoro/Voronoi3D.hpp
    fi
], [mpi_enabled=$mpi_enabled hdf5_enabled=$hdf5_enabled vtk_enabled=$vtk_enabled])

AC_LANG_POP()

AC_CONFIG_FILES([
    Makefile
    src/include/Makefile
    src/Makefile
])

AC_OUTPUT
