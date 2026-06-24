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
