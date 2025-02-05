AC_DEFUN([WITH_BOOST_CONFIGURE],
[
  DEFAULT_PATH=/usr
  DEFAULT_LIB_PATH=/usr/lib
  
  AC_ARG_WITH([boost],
              [AS_HELP_STRING([--with-boost(=DIR)],
              [The prefix where Boost is installed (default is the environment variable 'BOOST_DIR'), adding DIR/include])],
              [if test "$withval" = "yes"; then
                with_boost="no"
              else
                with_boost=$withval
              fi],
              [with_boost="no"]
             )

  AC_ARG_WITH([boost-libdir],
              [AS_HELP_STRING([--with-boost-libdir(=DIR)],
              ["Search path for boost libraries"])],
              [if test "$withval" = "yes"; then
                with_boost_libdir="no"
              else
                with_boost_libdir=$withval
              fi],
              [with_boost_libdir="no"]
             )

  AC_ARG_WITH([boost-include],
              [AS_HELP_STRING([--with-boost-include(=DIR)],
              ["Search path for boost header files"])],
              [if test "$withval" = "yes"; then
                with_boost_include="no"
              else
                with_boost_include=$withval
              fi],
              [with_boost_include="no"]
             )
])

AC_DEFUN([BOOST_CHECK_VERSION],
[
    major_required=`echo $1 | sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\1/'`
    minor_required=`echo $1 | sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\2/'`
    patch_required=`echo $1 | sed 's/\([[0-9]]*\).\([[0-9]]*\).\([[0-9]]*\)/\3/'`

    if test "$with_boost_include" = "no"; then
      AC_MSG_CHECKING([Boost's path])
      if test "$with_boost" = "no"; then
        with_boost=`echo $BOOST_DIR`
        if test -z "$with_boost"; then
          with_boost=`echo $BOOST_ROOT`
        fi
        if test -z "$with_boost"; then
          AC_MSG_ERROR([Boost not found])
        fi
      fi
      AC_MSG_RESULT($with_boost)
    fi

    if test "$with_boost_libdir" = "no"; then
      # Determine by default as the boost dir and then 'lib'
      if test -d "$with_boost"; then
        with_boost_libdir=`readlink -m $with_boost/lib`
      else
        AC_MSG_ERROR([Boost not found in $with_boost])
      fi
    fi

    if test "$with_boost_include" = "no"; then
      # Determine by default as the boost dir and then 'lib'
      if test -d "$with_boost"; then
        with_boost_include=`readlink -m $with_boost/include`
      else
        AC_MSG_ERROR([Boost not found in $with_boost])
      fi
    fi

    
    BOOST_CXXFLAGS="$CXXFLAGS -I $with_boost_include"
    BOOST_LDFLAGS="$LDFLAGS -L$with_boost_libdir -lboost_system"

    CXXFLAGS_BEFORE=$CXXFLAGS
    LDFLAGS_BEFORE=$LDFLAGS
    CXXFLAGS="$CXXFLAGS $BOOST_CXXFLAGS"
    LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
      
    # Check for the boost/version.hpp header file
    AC_CHECK_HEADER([boost/version.hpp], [have_boost=yes], [have_boost=no])
    
    AC_MSG_CHECKING([if boost's version is >= $1])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
                        [
                          #include <stdio.h>
                          #include <boost/version.hpp>
                        ],
                        [
                          constexpr size_t patch = BOOST_VERSION % 100;
                          constexpr size_t minor = (BOOST_VERSION / 100) % 100;
                          constexpr size_t major = BOOST_VERSION / 100000;
                          static_assert((major > $major_required) or 
                                        ((major == $major_required) and (minor > $minor_required)) or
                                        ((major == $major_required) and (minor == $minor_required) and (patch >= $patch_required)),
                                       "Boost version is too old");
                        ])
                      ], [boost_version_ok="yes"], [boost_version_ok="no"])

    if test "$boost_version_ok" = "yes"; then
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no])
      AC_MSG_ERROR([Boost version $1 or later required, but not found.])
    fi
])

AC_DEFUN([BOOST_INCLUDE_AND_LINK],
[
  CXXFLAGS=$BOOST_CXXFLAGS
  LDFLAGS=$BOOST_LDFLAGS
])

AC_DEFUN([BOOST_HANDLE],
[
  # first argument - version
  # second argument - what to do if success
  # third argument - what to do if failed

  WITH_BOOST_CONFIGURE
  BOOST_CHECK_VERSION([$1])
  if test "$boost_version_ok" = "yes" && test "$have_boost" = "yes"; then
    BOOST_INCLUDE_AND_LINK
    $2
  else
    $3
  fi
])