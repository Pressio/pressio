
# set where all the mains for google tests are
set(GTESTMAINSDIR "${CMAKE_CURRENT_SOURCE_DIR}/gtest_mains")
message(">> Google tests main files are in: " ${GTESTMAINSDIR})

# we need Gtest for unit tests
find_package(GTest REQUIRED)
if (NOT GTEST_FOUND)
message(FATAL_ERROR
	"Trying to build unit_tests but I cannot find GTest.
  	Please reconfigure using -DGTEST_ROOT=<path-to-gtest-install>.
  	The GTEST_ROOT should have a lib and include subdirs.")
endif()

include(macrosForCreatingUnitTests)
include(macrosForCreatingTests)

set(TESTING_LEVEL functional_small)

add_subdirectory(mpl)
add_subdirectory(utils)
add_subdirectory(type_traits)
add_subdirectory(expressions)
add_subdirectory(ops)
add_subdirectory(qr)
add_subdirectory(solvers_linear)
add_subdirectory(solvers_nonlinear)
add_subdirectory(ode)
