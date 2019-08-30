# define your minimum CMake version.
# This may be newer than the min required by TriBITS.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11 FATAL_ERROR)

# Make CMake set WIN32 with CYGWIN for older CMake versions.  CMake requires
# this to be in the top-level CMakeLists.txt file and not an include file :-(
#SET(CMAKE_LEGACY_CYGWIN_WIN32 1 CACHE BOOL "" FORCE)

set(pressio_ENABLE_CXX11 ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
MESSAGE("The c++ standard for PRESSIO is set for c++${CMAKE_CXX_STANDARD}" )

#
# A) Define your project name and set up major project options
#
# NOTE: Don't set options that would impact what packages get defined or
# enabled/disabled in this file as that would not impact other tools that
# don't process this file.
#

# Get PROJECT_NAME (must be in a file for other parts of system to read)
INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/ProjectName.cmake")

# CMake requires that you declare the CMake project in the
# top-level file and not in an include file :-(
PROJECT(${PROJECT_NAME} NONE)

#
# B) Pull in the TriBITS system and execute
#

SET(${PROJECT_NAME}_TRIBITS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/TriBITS/tribits" CACHE STRING "TriBITS base directory")

#SET(${PROJECT_NAME}_TRIBITS_DIR "${CMAKE_CURRENT_LIST_DIR}/../.."  CACHE  STRING
#  "TriBITS base directory (default assumes in TriBITS source tree)")
INCLUDE("${${PROJECT_NAME}_TRIBITS_DIR}/TriBITS.cmake")

# Set default location for header-only TPL to make easy to configure out of
# the TriBITS source tree.
SET(HeaderOnlyTpl_INCLUDE_DIRS
  "${${PROJECT_NAME}_TRIBITS_DIR}/examples/tpls/HeaderOnlyTpl"
  CACHE PATH "Default set by TriBITS/CMakeLists.txt" )

#------------------------------------------------------------
# set the directory where all the mains for google tests are
SET(GTESTMAINSDIR "${CMAKE_CURRENT_SOURCE_DIR}/gtest_mains")
MESSAGE("Google tests main files are in: " ${GTESTMAINSDIR})

#------------------------------------------------------------
# turn on cmake variables based on corresponding TPL
# these variables are defined inside the packages and
# are needed as preprocessor directives to turn on/off code

if (${CMAKE_CXX_STANDARD} EQUAL "14")
  MESSAGE("CXX_STANDARD=14, so HAVE_CPP14=ON")
  set(HAVE_CPP14 ON)
endif()

if( CMAKE_BUILD_TYPE STREQUAL "DEBUG" )
  message( "CMAKE_BUILD_TYPE==DEBUG so DEBUG_PRINT=ON" )
  set(DEBUG_PRINT ON)
endif()

if(TPL_ENABLE_GTEST)
  MESSAGE("TPL_ENABLE_GTEST=ON, so HAVE_GTEST=ON")
  set(HAVE_GTEST ON)
endif()

if(TPL_ENABLE_TRILINOS)
  MESSAGE("TPL_ENABLE_TRILINOS=ON, so HAVE_TRILINOS=ON")
  set(HAVE_TRILINOS ON)

  MESSAGE("TPL_ENABLE_TRILINOS=ON, so HAVE_KOKKOS=ON")
  set(HAVE_KOKKOS ON)

  if(DEBUG_PRINT)
    MESSAGE("TPL_ENABLE_TRILINOS=DEBUG_PRINT=ON so HAVE_TEUCHOS_TIMERS=ON")
    set(HAVE_TEUCHOS_TIMERS ON)
  endif()
endif()

if(TPL_ENABLE_KOKKOS)
  MESSAGE("TPL_ENABLE_KOKKOS=ON, so HAVE_KOKKOS=ON")
  set(HAVE_KOKKOS ON)
endif()

if(TPL_ENABLE_LAPACK)
  MESSAGE("TPL_ENABLE_LAPACK=ON, so HAVE_LAPACK=ON")
  set(HAVE_LAPACK ON)
endif()

if(TPL_ENABLE_PYBIND11)
  MESSAGE("TPL_ENABLE_PYBIND11=ON, so HAVE_PYBIND11=ON")
  set(HAVE_PYBIND11 ON)
  # if pybind11 is on, add include for python
  INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
endif()

if(TPL_ENABLE_BLAZE)
  MESSAGE("TPL_ENABLE_BLAZE=ON, so HAVE_BLAZE=ON")
  set(HAVE_BLAZE ON)
endif()

if(TPL_ENABLE_ARMADILLO)
  MESSAGE("TPL_ENABLE_ARMADILLO=ON, so HAVE_ARMADILLO=ON")
  set(HAVE_ARMADILLO ON)
endif()

if(DEBUG_PRINT)
  MESSAGE("DEBUG_PRINT is ON")
endif()


#------------------------------------------------------------------------
# this function creates and adds an individual unit-test executable
# currently it takes following arguments:
#
#  (1) test_name: self-explanatory
#
#  (2) mainfile: the name of the main file without extention.
#			  	 If file is ``test1.cc'', mainfile = test1
#
#  (3) ARGV2: the number of ranks, 1 or more
#
function(add_unit_test test_name mainfile )

  SET(NUMMPI ${ARGV2})
  if(${NUMMPI} EQUAL 1)
    SET(COMM1 serial mpi)
  elseif(${NUMMPI} GREATER 1)
    SET(COMM1 mpi)
  endif()
  # COMM serial => only run if pressio built without MPI
  # COMM mpi => only run if pressio built with MPI
  # COMM serial mpi => run in both cases,
  #		for MPI # of ranks set to tribits default 4

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    utest_${test_name}
    SOURCES ${test_name}.cc ${mainfile}.cc
    COMM ${COMM1}
    NUM_MPI_PROCS ${NUMMPI}
    PASS_REGULAR_EXPRESSION "PASSED"
    FAIL_REGULAR_EXPRESSION "FAILED"
  )
endfunction()

#------------------------------------------------------------------------
# this function creates and adds an individual SERIAL unit-test
# currently it takes following arguments:
#
#  (1) test_name: self-explanatory
#
function(add_serial_unit_test test_name)

  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    utest_${test_name}
    SOURCES ${test_name}.cc ${GTESTMAINSDIR}/gTestMain_serial.cc
    COMM serial mpi
    NUM_MPI_PROCS 1
    PASS_REGULAR_EXPRESSION "PASSED"
    FAIL_REGULAR_EXPRESSION "FAILED"
  )
endfunction()


#------------------------------------------------------------
# Do all of the processing for this Tribits project
TRIBITS_PROJECT()
