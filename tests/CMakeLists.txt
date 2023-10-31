include(CMakePrintHelpers)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(macrosForCreatingTests)
include(macrosForCreatingUnitTests)

# ---------------------------------
# 1. find or get gtest
# ---------------------------------
find_package(GTest)
if (GTEST_FOUND)
  add_library(gtest_main ALIAS GTest::gtest_main)
  cmake_print_variables(GTEST_INCLUDE_DIRS)
else()
  set(GTEST_VERSION "v1.14.0")
  set(BUILD_GMOCK   OFF)
  set(INSTALL_GTEST OFF)
  message(STATUS "GTest not found, fetching version ${GTEST_VERSION}")

  list(APPEND CMAKE_MESSAGE_INDENT "[GTest] ")
  include(FetchContent)
  FetchContent_Declare(
    googletest
    DOWNLOAD_EXTRACT_TIMESTAMP FALSE
    URL https://github.com/google/googletest/archive/refs/tags/${GTEST_VERSION}.tar.gz
    URL_HASH MD5=c8340a482851ef6a3fe618a082304cfc
  )
  FetchContent_MakeAvailable(googletest)
  list(POP_BACK CMAKE_MESSAGE_INDENT)
endif()

# ---------------------------------
# 2. where all mains are
# ---------------------------------
set(GTESTMAINSDIR "${CMAKE_CURRENT_SOURCE_DIR}/gtest_mains")
message("${Magenta}>> Google tests main files are in: ${GTESTMAINSDIR} ${ColourReset}")

# ---------------------------------
# 3. add test subdirectories
# ---------------------------------
if(PRESSIO_ENABLE_FUNCTIONAL_SMALL_TESTS OR PRESSIO_ENABLE_TESTS)
  add_subdirectory(functional_small)
endif()

if(PRESSIO_ENABLE_FUNCTIONAL_MEDIUM_TESTS OR PRESSIO_ENABLE_TESTS)
  add_subdirectory(functional_medium)
endif()

if(PRESSIO_ENABLE_FUNCTIONAL_LARGE_TESTS OR PRESSIO_ENABLE_TESTS)
  add_subdirectory(functional_large)
endif()
