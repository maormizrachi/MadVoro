function(add_madvoro_example target source)
  get_filename_component(_example_dir "${source}" DIRECTORY)

  add_executable(${target} ${source})
  target_link_libraries(${target} PRIVATE MadVoro::madvoro)
  target_include_directories(${target} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  set_target_properties(${target} PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/${_example_dir}"
  )
  add_custom_command(TARGET ${target} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "$<TARGET_FILE:${target}>"
      "${CMAKE_CURRENT_SOURCE_DIR}/${_example_dir}/"
    VERBATIM
  )

  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${_example_dir}/input")
    file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/${_example_dir}/input"
      DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/${_example_dir}"
    )
    add_custom_command(TARGET ${target} POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${CMAKE_CURRENT_SOURCE_DIR}/${_example_dir}/input"
        "$<TARGET_FILE_DIR:${target}>/input"
      VERBATIM
    )
  endif()
endfunction()

# Serial examples (no MPI required)
add_madvoro_example(example_faces_information faces_information/main.cpp)

add_madvoro_example(example_uniform_serial uniform_points/serial.cpp)

# Parallel examples (require MPI)
if(MADVORO_WITH_MPI)
  add_madvoro_example(example_uniform_parallel uniform_points/parallel.cpp)

  add_madvoro_example(example_pentagon pentagon/main.cpp)

  add_madvoro_example(example_pyramid pyramid/main.cpp)

  add_madvoro_example(example_fox fox/main.cpp)

  add_madvoro_example(example_parallel_multiple_builds parallel_multiple_builds/main.cpp)

  if(MADVORO_WITH_VTK)
    add_madvoro_example(example_parallel_output_vtk parallel_output_vtk/main.cpp)
  endif()
endif()

# ---------- CTest registration ----------
if(BUILD_TESTING)
  if(NOT MADVORO_WITH_MPI)
    add_test(NAME faces_information COMMAND example_faces_information)
    add_test(NAME uniform_serial    COMMAND example_uniform_serial)
    set_tests_properties(faces_information PROPERTIES
      WORKING_DIRECTORY "$<TARGET_FILE_DIR:example_faces_information>")
    set_tests_properties(uniform_serial PROPERTIES
      WORKING_DIRECTORY "$<TARGET_FILE_DIR:example_uniform_serial>")
  endif()

  if(MADVORO_WITH_MPI)
    set(_NP 4)
    add_test(NAME uniform_parallel
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_NP} $<TARGET_FILE:example_uniform_parallel>)
    set_tests_properties(uniform_parallel PROPERTIES
      WORKING_DIRECTORY "$<TARGET_FILE_DIR:example_uniform_parallel>")

    add_test(NAME pentagon
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_NP} $<TARGET_FILE:example_pentagon>)
    set_tests_properties(pentagon PROPERTIES
      WORKING_DIRECTORY "$<TARGET_FILE_DIR:example_pentagon>")

    add_test(NAME pyramid
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_NP} $<TARGET_FILE:example_pyramid>)
    set_tests_properties(pyramid PROPERTIES
      WORKING_DIRECTORY "$<TARGET_FILE_DIR:example_pyramid>")

    add_test(NAME fox
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_NP} $<TARGET_FILE:example_fox>)
    set_tests_properties(fox PROPERTIES
      WORKING_DIRECTORY "$<TARGET_FILE_DIR:example_fox>")
  endif()
endif()
