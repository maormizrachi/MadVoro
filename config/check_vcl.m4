AC_DEFUN([WITH_VCL_CONFIGURE],
[
  AC_ARG_WITH([vcl],
              [AS_HELP_STRING([--with-vcl(=DIR)],
              [The prefix where the VCL headers are installed])],
              [with_vcl=$withval],
              [with_vcl=no]
             )
])

AC_DEFUN([VCL_CHECK_VERSION],
[
    # first argument - version

    VCL_CXXFLAGS="-I $with_vcl"

    CXXFLAGS_BEFORE=$CXXFLAGS
    CXXFLAGS="$CXXFLAGS $VCL_CXXFLAGS"
    CPPFLAGS=$CXXFLAGS

    AC_CHECK_HEADER([vectorclass.h], [have_vcl=yes], [have_vcl=no])
    AC_MSG_CHECKING([if VCL's version is >= $1])

    AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
                        [
                          #include <vectorclass.h>
                        ],
                        [
                          static_assert(VECTORCLASS_H >= $1, "The current version of VCL (VECTORCLASS_H) is not compatible with required ($1)");
                        ])
                      ], [vcl_version_ok="yes"], [vcl_version_ok="no"])

    CXXFLAGS=$CXXFLAGS_BEFORE
    CPPFLAGS=$CXXFLAGS
    LDFLAGS=$LDFLAGS_BEFORE

    if test "$vcl_version_ok" = "yes"; then
      AC_MSG_RESULT([yes])
    else
      AC_MSG_RESULT([no])
    fi
])

AC_DEFUN([VCL_INCLUDE_AND_LINK],
[
  CXXFLAGS="$CXXFLAGS $VCL_CXXFLAGS"
  CPPFLAGS=$CXXFLAGS
])

AC_DEFUN([VCL_HANDLE],
[
  # first argument - version
  # second argument - what to do if success
  # third argument - what to do if failed

  WITH_VCL_CONFIGURE
  if test "$with_vcl" != "no"; then
    VCL_CHECK_VERSION([$1])
    if test "$vcl_version_ok" = "yes"; then
        VCL_INCLUDE_AND_LINK
        $2
    else
        $3
    fi
  else
    $3
  fi
])