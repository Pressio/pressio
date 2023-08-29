include(CMakePrintHelpers)

option(PRESSIO_ENABLE_CXX14 "Enable C++14" OFF)
option(PRESSIO_ENABLE_CXX17 "Enable C++17" OFF)
option(PRESSIO_ENABLE_CXX20 "Enable C++20" OFF)

option(PRESSIO_ENABLE_DEBUG_PRINT "Enable debug printing" OFF)
option(PRESSIO_ENABLE_TPL_EIGEN "Enable Eigen TPL" OFF)
option(PRESSIO_ENABLE_TPL_TRILINOS "Enable Trilinos TPL" OFF)
# option(PRESSIO_ENABLE_TEUCHOS_TIMERS "bla bla" OFF)
# option(PRESSIO_ENABLE_TPL_KOKKOS "bla bla" OFF)
# option(PRESSIO_ENABLE_TPL_MPI "bla bla" OFF)

if (CMAKE_CXX_STANDARD EQUAL 14)
  add_definitions(-DPRESSIO_ENABLE_CXX14)
elseif(CMAKE_CXX_STANDARD EQUAL 17)
  add_definitions(-DPRESSIO_ENABLE_CXX17)
elseif(CMAKE_CXX_STANDARD EQUAL 20)
  add_definitions(-DPRESSIO_ENABLE_CXX20)
endif()

if (PRESSIO_ENABLE_DEBUG_PRINT)
  add_definitions(-DPRESSIO_ENABLE_DEBUG_PRINT)
endif()

if(PRESSIO_ENABLE_TPL_EIGEN)
  if(NOT EIGEN_INCLUDE_DIR)
    message(FATAL_ERROR
      "I cannot find the Eigen headers. "
      "Please reconfigure with: -DEIGEN_INCLUDE_DIR=<full-path-to-headers>")
  endif()
  include_directories(${EIGEN_INCLUDE_DIR})
  add_definitions(-DPRESSIO_ENABLE_TPL_EIGEN)
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
  find_package(trilinos REQUIRED)
  add_definitions(-DPRESSIO_ENABLE_TPL_TRILINOS)

  # FINISH THIS
endif()

# if(PRESSIO_ENABLE_TPL_MPI)
#   find_package(MPI REQUIRED)
#   include_directories(${MPI_CXX_INCLUDE_DIRS})
# endif()

# if(PRESSIO_ENABLE_TPL_KOKKOS AND NOT PRESSIO_ENABLE_TPL_TRILINOS)
#   find_package(KokkosKernels REQUIRED)
#   set(KOKKOS_LIBS Kokkos::kokkoskernels)
#   #Kokkos::BLAS Kokkos::LAPACK Kokkos::kokkos Kokkos::kokkoskernels)
# endif()

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
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
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

# macro(add_utest_mpi TESTNAME TESTSRCS gMAIN nRANKS)
#   set(testNameFinal ${TESTNAME}_np${nRANKS})
#   add_executable(${testNameFinal} ${TESTSRCS} ${GTESTMAINSDIR}/${gMAIN}.cc)
#   target_link_libraries(${testNameFinal} ${MPI_CXX_LIBRARIES} gtest_main)
#   add_test(
#     NAME ${testNameFinal}
#     COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nRANKS}
#     ${MPIEXEC_PREFLAGS} ${testNameFinal} ${MPIEXEC_POSTFLAGS}
#     )
# endmacro()
#=====================================================================

# # macro to create and add an individual SERIAL test
# macro(add_serial_exe_and_test TESTNAME PKGNAME TESTSRCS REGEX)
#   set(testNameFinal ${PKGNAME}_${TESTNAME})
#   add_executable(${testNameFinal} ${TESTSRCS})
#   add_test(NAME ${testNameFinal} COMMAND ${testNameFinal})
#   set_tests_properties(
#     ${testNameFinal}
#     PROPERTIES PASS_REGULAR_EXPRESSION ${REGEX} FAIL_REGULAR_EXPRESSION "FAILED"
#   )
# endmacro()
# #=====================================================================

# macro(add_mpi_exe_and_test TESTNAME PKGNAME TESTSRCS nRANKS REGEXOK)
#   set(testNameFinal ${PKGNAME}_${TESTNAME}_np${nRANKS})
#   add_executable(${testNameFinal} ${TESTSRCS})
#   target_link_libraries(${testNameFinal} ${MPI_CXX_LIBRARIES})
#   add_test(
#     NAME ${testNameFinal}
#     COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nRANKS}
#     ${MPIEXEC_PREFLAGS} ${testNameFinal} ${MPIEXEC_POSTFLAGS}
#     )
#   set_tests_properties(
#     ${testNameFinal}
#     PROPERTIES PASS_REGULAR_EXPRESSION ${REGEXOK} FAIL_REGULAR_EXPRESSION "FAILED"
#   )
# endmacro()
# #=====================================================================

if(PRESSIO_ENABLE_FUNCTIONAL_SMALL_TESTS OR PRESSIO_ENABLE_TESTS)
  add_subdirectory(functional_small)
endif()

if(PRESSIO_ENABLE_FUNCTIONAL_MEDIUM_TESTS OR PRESSIO_ENABLE_TESTS)
  add_subdirectory(functional_medium)
endif()

if(PRESSIO_ENABLE_FUNCTIONAL_LARGE_TESTS OR PRESSIO_ENABLE_TESTS)
  add_subdirectory(functional_large)
endif()
