
# set where all the mains for google tests are
set(GTESTMAINSDIR "${CMAKE_CURRENT_SOURCE_DIR}/gtest_mains")
message("Google tests main files are in: " ${GTESTMAINSDIR})

# include the package directory where we have all the global headers
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../packages)

# we need Gtest for unit tests
find_package(GTest REQUIRED)
if (NOT GTEST_FOUND)
message(FATAL_ERROR
	"Trying to build unit_tests but I cannot find GTest.
  	Please reconfigure using -DGTEST_ROOT=<path-to-gtest-install>.
  	The GTEST_ROOT should have a lib and include subdirs.")
endif()

# load macros that help creating unit tests
include(macrosForCreatingUnitTests)

add_subdirectory(mpl)
add_subdirectory(utils)
add_subdirectory(containers)
add_subdirectory(ops)
add_subdirectory(qr)
add_subdirectory(optimizers)
add_subdirectory(solvers)
add_subdirectory(ode)
add_subdirectory(rom)
