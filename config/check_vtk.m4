AC_DEFUN([WITH_VTK_CONFIGURE],
[
  AC_ARG_WITH([vtk],
              [AS_HELP_STRING([--with-vtk(=DIR)],
              [The prefix where VTK is installed, adding DIR/include and DIR/lib])],
              [if test "$withval" = "yes"; then
                with_vtk=yes
              else
                with_vtk=$withval
              fi],
              [with_vtk=no]
             )

  AC_ARG_WITH([vtk-libdir],
              [AS_HELP_STRING([--with-vtk-libdir(=DIR)],
              ["Search path for VTK libraries"])],
              [if test "$withval" = "yes"; then
                with_vtk_libdir=no
              else
                with_vtk_libdir=$withval
              fi],
              [with_vtk_libdir=no]
             )

  AC_ARG_WITH([vtk-include],
              [AS_HELP_STRING([--with-vtk-libdir(=DIR)],
              ["Search path for VTK header files"])],
              [if test "$withval" = "yes"; then
                with_vtk_include=no
              else
                with_vtk_include=$withval
              fi],
              [with_vtk_include=no]
             )
])

AC_DEFUN([FIND_INSTALLED_VTK], [
    AC_MSG_CHECKING([Installed VTK path])
    commands=`compgen -c | grep "^vtk"`

    if test "$commands" = ""; then
      AC_MSG_ERROR([VTK could not be found])
    fi

    first_command_name=`compgen -c | grep "^vtk" | head -n 1`
    first_command=`which $first_command_name`
    first_dir=`dirname $first_command`

    # Flag to check if all directories match
    # Loop through the remaining commands and check their directories
    for cmd in $commands; do
        cmd_path=`which $cmd`
        cmd_dir=`dirname $cmd_path`
        if ! test "$cmd_dir" = "$first_dir"; then
            AC_MSG_ERROR([Two contradicting paths found: $first_dir (of command found at $first_command) and $cmd_dir (of command found at $command)])
            break
        fi
    done
    
    installed_directory=`realpath "$first_dir/.."`
    AC_MSG_RESULT([$installed_directory])
    with_vtk=$installed_directory
])

AC_DEFUN([VTK_CHECK_VERSION],
[
    # first argument - version, second argument - with mpi?

    if test "$with_vtk" = "yes"; then
      # '--with-vtk' but no path was given. Find this path
      FIND_INSTALLED_VTK
    fi

    if test "$with_vtk_include" = "no"; then
      with_vtk_include=`realpath $with_vtk/include`
    fi

    major_required=`echo $1 | sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\1/'`
    minor_required=`echo $1 | sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\2/'`
    build_required=`echo $1 | sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`
    version_dir=`ls $with_vtk_include`
    vtk_suffix="-`echo $version_dir | cut -d "-" -f 2`"
    vtk_include_files_dir=$with_vtk_include/$version_dir

    AC_MSG_CHECKING([VTK's header files location])
    AC_MSG_RESULT([$vtk_include_files_dir])

    AC_MSG_CHECKING([if VTK's version is >= $1])

    VTK_LIBS="-lvtkCommonCore$vtk_suffix -lvtkCommonColor$vtk_suffix -lvtkFiltersGeneral$vtk_suffix -lvtkFiltersSources$vtk_suffix -lvtkIOXML$vtk_suffix -lvtkInteractionStyle$vtk_suffix"
    if test "$2" = "yes"; then
        VTK_LIBS="$VTK_LIBS -lvtkIOParallelXML$vtk_suffix -lvtkParallelMPI$vtk_suffix"
    fi
    
    VTK_CXXFLAGS="-I $vtk_include_files_dir"
    VTK_LDFLAGS="-L$with_vtk/lib -L$with_vtk/lib64"
    if test "$with_vtk_libdir" != "no"; then
        VTK_LDFLAGS="$VTK_LDFLAGS -L$with_vtk_libdir"
    fi
    VTK_LDFLAGS="$VTK_LDFLAGS $VTK_LIBS"

    CXXFLAGS_BEFORE=$CXXFLAGS
    CXXFLAGS="$CXXFLAGS $VTK_CXXFLAGS"
    CPPFLAGS=$CXXFLAGS
    LDFLAGS_BEFORE=$LDFLAGS
    LDFLAGS="$LDFLAGS $VTK_LDFLAGS"
    
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
                        [
                          #include <stdio.h>
                          #include <vtkVersionMacros.h>

                          #define VTK_VER (VTK_MAJOR_VERSION.VTK_MINOR_VERSION.VTK_BUILD_VERSION) 
                        ],
                        [
                          static_assert((VTK_MAJOR_VERSION > $major_required)
                                        or ((VTK_MAJOR_VERSION == $major_required) and (VTK_MINOR_VERSION >= $minor_required))
                                        or ((VTK_MAJOR_VERSION == $major_required) and (VTK_MINOR_VERSION == $minor_required) and (VTK_BUILD_VERSION >= $build_required)),
                                        "VTK version (VTK_VER) is too old");
                        ])
                      ], [vtk_version_ok="yes"], [vtk_version_ok="no"])

    CXXFLAGS=$CXXFLAGS_BEFORE
    CPPFLAGS=$CXXFLAGS
    LDFLAGS=$LDFLAGS_BEFORE

    if test "$vtk_version_ok" = "yes"; then
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([VTK version $1 or later required, but not found. Change version, or compile without VTK support.])
    fi
])

AC_DEFUN([VTK_INCLUDE_AND_LINK],
[
  CXXFLAGS="$CXXFLAGS $VTK_CXXFLAGS"
  CPPFLAGS=$CXXFLAGS
  LDFLAGS="$LDFLAGS $VTK_LDFLAGS"
])

AC_DEFUN([VTK_HANDLE],
[
  # first argument - version
  # second argument - with mpi
  # third argument - what to do if success
  # fourth argument - what to do if failed

  WITH_VTK_CONFIGURE
  if test "$with_vtk" != "no" || test "$with_vtk_libdir" != "no"; then
    VTK_CHECK_VERSION([$1], [$2])
    if test "$vtk_version_ok" = "yes"; then
        VTK_INCLUDE_AND_LINK
        $3
    else
        $4
    fi
  else
    $4
  fi
])