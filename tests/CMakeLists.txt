
include_directories(
  ${CMAKE_BINARY_DIR} # needed for the pressio_cmake_config.h
  ${CMAKE_SOURCE_DIR}/include # include all the pressio headers
  ${CMAKE_CURRENT_SOURCE_DIR}/testing_apps
  )

# load macros that help creating tests
include(macrosForCreatingTests)

add_subdirectory(optimizers)
add_subdirectory(solvers_linear)
add_subdirectory(solvers_nonlinear)

#add_subdirectory(ode)
#add_subdirectory(rom)
#add_subdirectory(apps)
