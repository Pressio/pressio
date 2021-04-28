
include_directories(
  ${CMAKE_BINARY_DIR} # needed for the pressio_cmake_config.h
  ${CMAKE_SOURCE_DIR}/include # include all the global headers
  )

# load macros that help creating tests
include(macrosForCreatingTests)

add_subdirectory(containers)
add_subdirectory(apps)
add_subdirectory(solvers)
add_subdirectory(ode)
add_subdirectory(optimizers)
add_subdirectory(rom)
