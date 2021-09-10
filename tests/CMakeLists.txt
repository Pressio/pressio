
include_directories(${CMAKE_BINARY_DIR} ${CMAKE_SOURCE_DIR}/include)

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

if(PRESSIO_ENABLE_FUNCTIONAL_SMALL_TESTS OR PRESSIO_ENABLE_TESTS)
  add_subdirectory(functional_small)
endif()

if(PRESSIO_ENABLE_FUNCTIONAL_MEDIUM_TESTS OR PRESSIO_ENABLE_TESTS)
  add_subdirectory(functional_medium)
endif()

if(PRESSIO_ENABLE_FUNCTIONAL_LARGE_TESTS OR PRESSIO_ENABLE_TESTS)
  add_subdirectory(functional_large)
endif()
