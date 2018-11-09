################################################################################
#                                                                              #
#                        TribitsExampleProjectMeta                             #
#                                                                              #
################################################################################

# To be safe, define your minimum CMake version.  This may be newer than the
# min required by TriBITS.
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.11 FATAL_ERROR)

# Make CMake set WIN32 with CYGWIN for older CMake versions.  CMake requires
# this to be in the top-level CMakeLists.txt file and not an include file :-(
SET(CMAKE_LEGACY_CYGWIN_WIN32 1 CACHE BOOL "" FORCE)

set(rompp_ENABLE_CXX11 ON)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)
MESSAGE("The c++ standard is set for c++${CMAKE_CXX_STANDARD}" )

#
# A) Define your project name and set up major project options
#
# NOTE: Don't set options that would impact what packages get defined or
# enabled/disabled in this file as that would not impact other tools that
# don't process this file.
#

# Get PROJECT_NAME (must be in a file for other parts of system to read)
INCLUDE("${CMAKE_CURRENT_SOURCE_DIR}/ProjectName.cmake")

# CMake requires that you declare the CMake project in the top-level file and
# not in an include file :-(
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
# these variables are defined inside the core package and 
# are needed as preprocessor directives to turn on/off code
# if the TPL is not linked to

if(TPL_ENABLE_BLAZE)
  set(HAVE_BLAZE ON)
endif()
if(TPL_ENABLE_MPI)
  set(HAVE_MPI ON)
endif()
if(TPL_ENABLE_TRILINOS)
  set(HAVE_TRILINOS ON)
endif()
if(TPL_ENABLE_ARMADILLO)
  set(HAVE_ARMADILLO ON)
endif()

if (${CMAKE_CXX_STANDARD} EQUAL "14")
  set(HAVE_CPP14 ON)
endif()

#------------------------------------------------------------
# Do all of the processing for this Tribits project
TRIBITS_PROJECT()
