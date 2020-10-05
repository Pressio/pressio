
# include the package directory where we have all the global headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../packages)

# load macros that help creating tests
include(macrosForCreatingTests)

add_subdirectory(containers)
add_subdirectory(apps)
add_subdirectory(solvers)
add_subdirectory(ode)
add_subdirectory(optimizers)
add_subdirectory(rom)
