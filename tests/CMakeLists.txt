
include_directories(
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/include
)

# load macros that help creating tests
include(macrosForCreatingTests)

add_subdirectory(solvers_nonlinear)

#add_subdirectory(optimizers)
#add_subdirectory(ode)
#add_subdirectory(rom)
#add_subdirectory(apps)
