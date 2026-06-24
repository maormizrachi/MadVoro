# Serial examples (no MPI required)
add_executable(example_faces_information faces_information/main.cpp)
target_link_libraries(example_faces_information PRIVATE MadVoro::madvoro)
target_include_directories(example_faces_information PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(example_uniform_serial uniform_points/serial.cpp)
target_link_libraries(example_uniform_serial PRIVATE MadVoro::madvoro)
target_include_directories(example_uniform_serial PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Parallel examples (require MPI)
if(MADVORO_WITH_MPI)
  add_executable(example_uniform_parallel uniform_points/parallel.cpp)
  target_link_libraries(example_uniform_parallel PRIVATE MadVoro::madvoro)
  target_include_directories(example_uniform_parallel PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

  add_executable(example_pentagon pentagon/main.cpp)
  target_link_libraries(example_pentagon PRIVATE MadVoro::madvoro)
  target_include_directories(example_pentagon PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

  add_executable(example_pyramid pyramid/main.cpp)
  target_link_libraries(example_pyramid PRIVATE MadVoro::madvoro)
  target_include_directories(example_pyramid PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

  add_executable(example_fox fox/main.cpp)
  target_link_libraries(example_fox PRIVATE MadVoro::madvoro)
  target_include_directories(example_fox PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

  add_executable(example_parallel_multiple_builds parallel_multiple_builds/main.cpp)
  target_link_libraries(example_parallel_multiple_builds PRIVATE MadVoro::madvoro)
  target_include_directories(example_parallel_multiple_builds PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

  if(MADVORO_WITH_VTK)
    add_executable(example_parallel_output_vtk parallel_output_vtk/main.cpp)
    target_link_libraries(example_parallel_output_vtk PRIVATE MadVoro::madvoro)
    target_include_directories(example_parallel_output_vtk PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  endif()
endif()

# ---------- CTest registration ----------
if(BUILD_TESTING)
  if(NOT MADVORO_WITH_MPI)
    add_test(NAME faces_information COMMAND example_faces_information)
    add_test(NAME uniform_serial    COMMAND example_uniform_serial)
  endif()

  if(MADVORO_WITH_MPI)
    set(_NP 4)
    add_test(NAME uniform_parallel
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_NP} $<TARGET_FILE:example_uniform_parallel>)

    add_test(NAME pentagon
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_NP} $<TARGET_FILE:example_pentagon>)

    add_test(NAME pyramid
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_NP} $<TARGET_FILE:example_pyramid>)

    add_test(NAME fox
      COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${_NP} $<TARGET_FILE:example_fox>)
    set_tests_properties(fox PROPERTIES
      WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/fox")
  endif()
endif()
