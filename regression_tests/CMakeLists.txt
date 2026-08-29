# Build one executable for every MadVoro regression case.
#
# REGRESSION_INFO is the source of truth for THUNDER discovery.  Each case
# directory must contain REGRESSION_INFO and a main source file.

set(_madvoro_regression_info_files "")
file(GLOB _madvoro_regression_info_files CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/cases/*/REGRESSION_INFO")
list(SORT _madvoro_regression_info_files)

set(_madvoro_regression_include_dirs
  "${CMAKE_CURRENT_SOURCE_DIR}/.."
  "${CMAKE_CURRENT_SOURCE_DIR}/lib"
  "${CMAKE_CURRENT_SOURCE_DIR}/../examples")

foreach(_madvoro_info_file IN LISTS _madvoro_regression_info_files)
  get_filename_component(_madvoro_case_dir "${_madvoro_info_file}" DIRECTORY)
  get_filename_component(_madvoro_case_name "${_madvoro_case_dir}" NAME)

  file(STRINGS "${_madvoro_info_file}" _madvoro_tags_lines
    REGEX "^[ \\t]*TAGS[ \\t]*=")
  string(TOLOWER "${_madvoro_tags_lines}" _madvoro_tags)
  if(_madvoro_tags MATCHES "mpi")
    set(_madvoro_case_requires_mpi ON)
  else()
    set(_madvoro_case_requires_mpi OFF)
  endif()

  if(_madvoro_case_requires_mpi AND NOT MADVORO_WITH_MPI)
    message(STATUS
      "Skipping ${_madvoro_case_name}: REGRESSION_INFO marks it as MPI "
      "but this build has MPI disabled")
    continue()
  endif()

  file(GLOB _madvoro_case_sources CONFIGURE_DEPENDS
    "${_madvoro_case_dir}/*.c"
    "${_madvoro_case_dir}/*.cc"
    "${_madvoro_case_dir}/*.cpp"
    "${_madvoro_case_dir}/*.cxx")
  list(SORT _madvoro_case_sources)

  set(_madvoro_case_has_main OFF)
  foreach(_madvoro_case_source IN LISTS _madvoro_case_sources)
    get_filename_component(_madvoro_case_source_name
      "${_madvoro_case_source}" NAME)
    if(_madvoro_case_source_name MATCHES "^main\\.(c|cc|cpp|cxx)$")
      set(_madvoro_case_has_main ON)
    endif()
  endforeach()

  if(NOT _madvoro_case_has_main)
    message(WARNING
      "Skipping ${_madvoro_case_name}: ${_madvoro_info_file} exists, but "
      "the directory has no main source file")
    continue()
  endif()

  if(TARGET "${_madvoro_case_name}")
    message(FATAL_ERROR
      "Cannot create MadVoro regression target ${_madvoro_case_name}: "
      "the target name is already in use")
  endif()

  add_executable("${_madvoro_case_name}" ${_madvoro_case_sources})
  target_link_libraries("${_madvoro_case_name}" PRIVATE MadVoro::madvoro)
  target_include_directories("${_madvoro_case_name}" PRIVATE
    ${_madvoro_regression_include_dirs})
  target_compile_definitions("${_madvoro_case_name}" PRIVATE
    "MADVORO_DATA_DIR=\"${_madvoro_case_dir}\"")
  set_target_properties("${_madvoro_case_name}" PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${_madvoro_case_dir}")

  if(NOT _madvoro_case_requires_mpi AND MADVORO_WITH_MPI)
    target_compile_options("${_madvoro_case_name}" PRIVATE
      -UMADVORO_WITH_MPI -USPATIAL_DS_WITH_MPI -U__WITH_MPI)
  endif()

  message(STATUS "MadVoro regression case: ${_madvoro_case_name}")
endforeach()
