include(CMakePrintHelpers)

# ---------------------------------
# 1. find or get gtest
# ---------------------------------
# adapted from Kokkos
find_package(GTest)
if (GTEST_FOUND)
  add_library(gtest_main ALIAS GTest::gtest_main)
  cmake_print_variables(GTEST_INCLUDE_DIRS)
else()
  set(GTEST_VERSION "v1.14.0")
  message(STATUS "GTest not found, fetching version ${GTEST_VERSION}")

  list(APPEND CMAKE_MESSAGE_INDENT "[GTest] ")
  include(FetchContent)
  FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
    URL_HASH MD5=c8340a482851ef6a3fe618a082304cfc
    DOWNLOAD_EXTRACT_TIMESTAMP FALSE
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
# 3. macros for createing unit tests
# ---------------------------------
macro(add_serial_utest TESTNAME)
  set(testNameFinal ${TESTNAME})
  add_executable(${testNameFinal} ${ARGN} ${GTESTMAINSDIR}/gTestMain_serial.cc)
  target_link_libraries(${testNameFinal} pressio gtest_main)
  add_test(NAME ${testNameFinal} COMMAND ${testNameFinal})
endmacro()

macro(add_serial_utest_kokkos TESTNAME TESTSRCS)
  set(testNameFinal ${TESTNAME})
  add_executable(${testNameFinal} ${TESTSRCS} ${GTESTMAINSDIR}/gTestMain_kokkos.cc)
  target_link_libraries(${testNameFinal} ${KOKKOS_LIBS} pressio gtest_main)
  add_test(NAME ${testNameFinal} COMMAND ${testNameFinal})
endmacro()

macro(add_utest_mpi TESTNAME TESTSRCS gMAIN nRANKS)
  set(testNameFinal ${TESTNAME}_np${nRANKS})
  add_executable(${testNameFinal} ${TESTSRCS} ${GTESTMAINSDIR}/${gMAIN}.cc)
  target_link_libraries(${testNameFinal} ${MPI_CXX_LIBRARIES} pressio gtest_main)
  add_test(
    NAME ${testNameFinal}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nRANKS}
    ${MPIEXEC_PREFLAGS} ${testNameFinal} ${MPIEXEC_POSTFLAGS}
    )
endmacro()
#=====================================================================

if(PRESSIO_ENABLE_FUNCTIONAL_SMALL_TESTS OR PRESSIO_ENABLE_TESTS)
  add_subdirectory(functional_small)
endif()

if(PRESSIO_ENABLE_FUNCTIONAL_MEDIUM_TESTS OR PRESSIO_ENABLE_TESTS)
  add_subdirectory(functional_medium)
endif()

if(PRESSIO_ENABLE_FUNCTIONAL_LARGE_TESTS OR PRESSIO_ENABLE_TESTS)
  add_subdirectory(functional_large)
endif()
