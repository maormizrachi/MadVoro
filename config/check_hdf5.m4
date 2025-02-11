AC_DEFUN([WITH_HDF5_CONFIGURE],
[
  AC_ARG_WITH([hdf5],
              [AS_HELP_STRING([--with-hdf5(=DIR)],
              [The prefix where HDF5 With CPP support is installed (default (no DIR) is automatic finding), adding DIR/bin, DIR/include and DIR/lib])],
              [
                if test "$withval" = "yes"; then
                    with_hdf5=yes
                else
                    with_hdf5=$withval
                fi
              ],
              [with_hdf5=no]
             )
    

  AC_ARG_WITH([hdf5-libdir],
              [AS_HELP_STRING([--with-hdf5-libdir],
              ["Search path for HDF5 With CPP support libraries"])],
              [
                if test "$withval" = "yes"; then
                    with_hdf5_libdir=no
                else
                    with_hdf5_libdir=$withval
                fi
              ],
              [with_hdf5_libdir=no]
             )


  AC_ARG_WITH([hdf5-include],
              [AS_HELP_STRING([--with-hdf5-include],
              ["Search path for HDF5 With CPP support header files"])],
              [
                if test "$withval" = "yes"; then
                    with_hdf5_include=no
                else
                    with_hdf5_include=$withval
                fi
              ],
              [with_hdf5_include=no]
             )
    
])

AC_DEFUN([HDF5_CHECK_VERSION],
[
    # first argument - version, second argument - with mpi

    if test "$with_hdf5" = "yes"; then
      h5dump=`command -v h5dump`
      hdf5_bin_dir=`dirname $h5dump`
      with_hdf5=`realpath "$hdf5_bin_dir/.."`
    fi

    if test "$with_hdf5_libdir" = "no"; then 
      # did not supply any value with 'with-hdf5-libdir', get by path
      with_hdf5_libdir=`realpath "$with_hdf5/lib"`
    fi

    if test "$with_hdf5_include" = "no"; then 
      # did not supply any value with 'with-hdf5-include', get by path
      with_hdf5_include=`realpath "$with_hdf5/include"`
    fi

    AC_MSG_CHECKING([HDF5 path])
    AC_MSG_RESULT([$with_hdf5])

    major_required=`echo $1 | sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\1/'`
    minor_required=`echo $1 | sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\2/'`
    release_required=`echo $1 | sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`
    
    HDF5_LIBS="-lhdf5_cpp -lhdf5_hl_cpp"
    HDF5_CXXFLAGS="$CXXFLAGS -I $with_hdf5_include"
    HDF5_LDFLAGS="$LDFLAGS -L$with_hdf5_libdir"
    HDF5_LDFLAGS="$HDF5_LDFLAGS $HDF5_LIBS"

    AC_MSG_CHECKING([HDF5's header files location])
    AC_MSG_RESULT([$with_hdf5_include])
    AC_CHECK_HEADER([H5public.h], [have_h5_header=yes], [have_h5_header=no])

    CXXFLAGS_BEFORE=$CXXFLAGS
    LDFLAGS_BEFORE=$LDFLAGS
    CXXFLAGS="$CXXFLAGS $HDF5_CXXFLAGS"
    CPPFLAGS=$CXXFLAGS
    LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"

    AC_MSG_CHECKING([if HDF5's version is >= $1])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
                        [
                          #include <H5public.h>

                          #define HDF5_VER (H5_VERS_MAJOR.H5_VERS_MINOR.H5_VERS_RELEASE) 
                        ],
                        [
                          static_assert((H5_VERS_MAJOR > $major_required)
                                        or ((H5_VERS_MAJOR == $major_required) and (H5_VERS_MINOR >= $minor_required))
                                        or ((H5_VERS_MAJOR == $major_required) and (H5_VERS_MINOR == $minor_required) and (H5_VERS_RELEASE >= $release_required)),
                                        "HDF5 version (HDF5_VER) is too old");
                        ])
                      ], [hdf5_version_ok="yes"], [hdf5_version_ok="no"])

    CXXFLAGS="$CXXFLAGS_BEFORE"
    CPPFLAGS=$CXXFLAGS
    LDFLAGS="$LDFLAGS_BEFORE"

    if test "$hdf5_version_ok" = "yes"; then
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([HDF5 version $1 or later required, but not found. Change version, or compile without HDF5 support.])
    fi
])

AC_DEFUN([HDF5_CHECK_CPP],
[    
    CXXFLAGS_BEFORE=$CXXFLAGS
    LDFLAGS_BEFORE=$LDFLAGS
    CXXFLAGS="$CXXFLAGS $HDF5_CXXFLAGS"
    LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"

    CPPFLAGS_BEFORE=$CPPFLAGS
    CPPFLAGS=$CXXFLAGS

    # Check for the H5Cpp.h header file
    AC_CHECK_HEADER([H5Cpp.h], [have_hdf5_cpp=yes], [have_hdf5_cpp=no])
    
    # Check for the C++ library
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
                        [
                          #include <string>
                          #include <H5Cpp.h>

                          void test_H5std_string()
                          {
                            H5std_string s("test");
                          }
                        ],
                        [
                          test_H5std_string();
                          return 0;
                        ])
                      ], [hdf5_cpp_ok="yes"], [hdf5_cpp_ok="no"])

    CXXFLAGS=$CXXFLAGS_BEFORE
    CPPFLAGS=$CPPFLAGS_BEFORE
    LDFLAGS=$LDFLAGS_BEFORE
])

AC_DEFUN([HDF5_INCLUDE_AND_LINK],
[
  CXXFLAGS="$CXXFLAGS $HDF5_CXXFLAGS"
  CPPFLAGS=$CXXFLAGS
  LDFLAGS="$LDFLAGS $HDF5_LDFLAGS"
])

AC_DEFUN([HDF5_HANDLE],
[
  # first argument - version
  # second argument - what to do if success
  # third argument - what to do if failed

  WITH_HDF5_CONFIGURE
  if test "$with_hdf5" != "no" || test "$with_hdf5_libdir" != "no"; then
    HDF5_CHECK_VERSION([$1])
    if test "$hdf5_version_ok" = "yes"; then
        HDF5_CHECK_CPP
        if test "$hdf5_cpp_ok" = "yes"; then
            HDF5_INCLUDE_AND_LINK
            $2
        else
            $3
        fi
    else
        $3
    fi
  else
    $3
  fi
])