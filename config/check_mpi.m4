AC_DEFUN([WITH_MPI_CONFIGURE],
[
  AC_ARG_WITH([mpi],
              [AS_HELP_STRING([--with-mpi(=MPIHOME)], [Enable MPI support. mpicxx should be in MPIHOME/bin, include files should be in MPIHOME/include and libraries in MPIHOME/lib])],
              [if test "$withval" = "yes"; then
                with_mpi=yes
              else
                with_mpi=$withval
              fi],
              [with_mpi=no])

  # Define variables for MPI compilers, allowing user to specify them
  AC_ARG_VAR([MPICXX], [C++ compiler for MPI (e.g., mpicxx, mpiicpc)])
])

AC_DEFUN([MPI_CHECK_VERSION],
[
    if test "$MPICXX" = ""; then
      # Only if did not supply any MPICXX value

      if test "$with_mpi" = "yes"; then
        # did not supply any value with 'with-mpi', get by path

        MPICXX=`command -v mpicxx`
        mpi_bin_dir=`dirname $MPICXX`
        with_mpi="$mpi_bin_dir/.."

      else
        # did not supply value for MPICXX, but did supply 'with-mpi'

        MPICXX="`realpath $with_mpi`/bin/mpicxx"
      fi
    fi

    # Check MPICXX
    AC_MSG_CHECKING(MPI compiler path)
    AC_MSG_RESULT($MPICXX)
    
    if ! test -n "$MPICXX"; then
        AC_MSG_ERROR([MPI compilers not found])
    fi

    major_required=`echo $1 | sed 's/\([[0-9]]*\).\([[0-9]]*\)/\1/'`
    minor_required=`echo $1 | sed 's/\([[0-9]]*\).\([[0-9]]*\)/\2/'`

    CXX_BEFORE=$CXX
    CXXFLAGS_BEFORE=$CXXFLAGS
    LDFLAGS_BEFORE=$LDFLAGS

    MPI_CXXFLAGS=""
    MPI_LDFLAGS=""

    CXX=$MPICXX
    CXXFLAGS="$CXXFLAGS $MPI_CXXFLAGS"
    CPPFLAGS=$CXXFLAGS
    LDFLAGS="$LDFLAGS $MPI_LDFLAGS"

    AC_CHECK_HEADER([mpi.h], [have_mpi_header=yes], [have_mpi_header=no])

    AC_MSG_CHECKING([if mpi's version is >= $1])

    AC_RUN_IFELSE([AC_LANG_PROGRAM(
                    [
                      #include <iostream>
                      #include <mpi.h>
                    ],
                    [
                      MPI_Init(NULL, NULL);
                      int major, minor;
                      MPI_Get_version(&major, &minor);
                      int ok = 1;
                      if((major > $major_required) or ((major == $major_required) and (minor >= $minor_required)))
                      {
                          ok = 1;
                      }
                      MPI_Finalize();
                      return (1 - ok);
                    ])
                  ], [mpi_version_ok="yes"], [mpi_version_ok="no"])

    CXXFLAGS=$CXXFLAGS_BEFORE
    CPPFLAGS=$CXXFLAGS
    LDFLAGS=$LDFLAGS_BEFORE
    CXX=$CXX_BEFORE

    if test "$mpi_version_ok" = "yes"; then
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([MPI supporting standard $1 or later is required, but not found. Change version, or compile without MPI support.])
    fi
])

AC_DEFUN([MPI_INCLUDE_AND_LINK],
[
  CXX=$MPICXX
  CXXFLAGS="$CXXFLAGS $MPI_CXXFLAGS"
  CPPFLAGS=$CXXFLAGS
  LDFLAGS="$LDFLAGS $MPI_LDFLAGS"
])

AC_DEFUN([MPI_HANDLE],
[
  # first argument - version
  # second argument - what to do if success
  # third argument - what to do if failed

  WITH_MPI_CONFIGURE
  if test "$with_mpi" != "no"; then
    MPI_CHECK_VERSION([$1])
    if test "$mpi_version_ok" != "no"; then
      MPI_INCLUDE_AND_LINK
      $2
    else
      $3
    fi
  fi
])