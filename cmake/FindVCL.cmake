# FindVCL.cmake - Find the Vector Class Library (header-only)
#
# Input variables:
#   VCL_DIR  - hint directory containing vectorclass.h
#
# Output variables:
#   VCL_FOUND        - TRUE if found
#   VCL_INCLUDE_DIRS - include directory

find_path(VCL_INCLUDE_DIR
  NAMES vectorclass.h
  HINTS
    ${VCL_DIR}
    $ENV{VCL_DIR}
    $ENV{VCL_ROOT}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(VCL
  REQUIRED_VARS VCL_INCLUDE_DIR
)

if(VCL_FOUND)
  set(VCL_INCLUDE_DIRS ${VCL_INCLUDE_DIR})
endif()

mark_as_advanced(VCL_INCLUDE_DIR)
