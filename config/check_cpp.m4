# Check for C++ standard flag support
AC_MSG_CHECKING([for C++17 support])

CXXFLAGS="$CXXFLAGS -std=c++17"

AC_COMPILE_IFELSE([AC_LANG_PROGRAM([], [[
#if __cplusplus < 201703L
#error "C++17 is required"
#endif
]])], [
  AC_MSG_RESULT([yes])
], [
  AC_MSG_ERROR([C++17 support is required. Use a compiler that supports C++17 or later.])
])

