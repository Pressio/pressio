
cmake_minimum_required(VERSION 3.11.0 FATAL_ERROR)
project(pressio C CXX)

# setting where other cmake commands are
set(CMAKE_MODULE_PATH
  ${PROJECT_SOURCE_DIR}/cmake;${PROJECT_SOURCE_DIR}/cmake/tpls)

# include commands for coloring prints
include(colors)

#=====================================================================
# versioning
#=====================================================================
set(Pressio_VERSION_MAJOR 0)
set(Pressio_VERSION_MINOR 8)
set(Pressio_VERSION_PATCH 0)
set(Pressio_VERSION
  "${Pressio_VERSION_MAJOR}.${Pressio_VERSION_MINOR}.${Pressio_VERSION_PATCH}")
math(EXPR PRESSIO_VERSION
  "${Pressio_VERSION_MAJOR} * 10000 + ${Pressio_VERSION_MINOR} * 100 + ${Pressio_VERSION_PATCH}")

message("${Yellow}----------------------------------${ColourReset}")
message("${Yellow} Pressio version = ${Pressio_VERSION} ${ColourReset}")
message("${Yellow}----------------------------------${ColourReset}")

#=====================================================================
# guard here for specific cases
#=====================================================================
# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed.
    Please make a new directory (called a build directory) and run CMake from there.
    You may need to remove CMakeCache.txt.")
endif()

# default to relase if build type is empty
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if( NOT cmake_build_type_tolower STREQUAL "debug"
    AND NOT cmake_build_type_tolower STREQUAL "release")
  message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\".
    Allowed values are Debug, Release (case-insensitive).")
endif()

#=====================================================================
# c++ standard
#=====================================================================
# commands to test if compiler supports standard
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORT_CPP11)
if(NOT COMPILER_SUPPORT_CPP11)
  message(FATAL_ERROR "Compiler does not support -std=c++11. This is required.")
endif()
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#=====================================================================
# other stuff
#=====================================================================
include(GNUInstallDirs)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_subdirectory(packages)

#=====================================================================
# note that actually linking the TPLs libraries is only needed
# when building tests or when the user builds their app.
# For just installing pressio, we only need to set the cmake vars.
#=====================================================================
message("")
message("${Magenta}Processing what is being enabled/disabled:${ColourReset}")

if( cmake_build_type_tolower STREQUAL "debug" AND NOT PRESSIO_ENABLE_DEBUG_PRINT)
  set(PRESSIO_ENABLE_DEBUG_PRINT ON)
  message( "> CMAKE_BUILD_TYPE==Debug, setting PRESSIO_ENABLE_DEBUG_PRINT=${PRESSIO_ENABLE_DEBUG_PRINT}" )
  message( "  To disable: reconfigure with -DPRESSIO_ENABLE_DEBUG_PRINT=OFF" )
  message("")
endif()
if(PRESSIO_ENABLE_DEBUG_PRINT)
  message( "> Enabling debug prints since PRESSIO_ENABLE_DEBUG_PRINT=ON" )
  message( "  To disable: reconfigure with -DPRESSIO_ENABLE_DEBUG_PRINT=OFF" )
  message("")
endif()

option(PRESSIO_ENABLE_UNIT_TESTS "Enabling pressio unit testing." OFF)
option(PRESSIO_ENABLE_TESTS "Enable ALL pressio tests: regression and unit tests." OFF)

option(PRESSIO_ENABLE_TPL_EIGEN		"Enable Eigen TPL"	ON)
option(PRESSIO_ENABLE_TPL_TRILINOS	"Enable Trilinos TPL"	OFF)
option(PRESSIO_ENABLE_TPL_BLAS		"Enable BLAS TPL"	OFF)
option(PRESSIO_ENABLE_TPL_LAPACK	"Enable LAPACK TPL"	OFF)
option(PRESSIO_ENABLE_TPL_KOKKOS	"Enable Kokkos TPL"	OFF)
option(PRESSIO_ENABLE_TPL_MPI		"Enable MPI"		OFF)
option(PRESSIO_ENABLE_TPL_PYBIND11	"Enable Pybind11 TPL"	OFF)

include(tplEigen)
include(tplPybind11)
include(tplMPI)
include(tplTrilinos)
include(tplBlas)
include(tplLapack)
include(tplKokkos)
message("")

#=====================================================================
# unit tests are enabled if
#   (a) PRESSIO_ENABLE_UNIT_TESTS = ON
#   (b) PRESSIO_ENABLE_TESTS = ON
#=====================================================================
if(PRESSIO_ENABLE_UNIT_TESTS OR PRESSIO_ENABLE_TESTS)
  enable_testing()
  add_subdirectory(unit_tests)
endif()

#=====================================================================
# basic regression tests are enabled if
#   PRESSIO_ENABLE_TESTS = ON
#=====================================================================
if(PRESSIO_ENABLE_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

#=====================================================================
# configure file with some cmakedefines that are used
# within the code in preprocessor directives statements
#=====================================================================
configure_file(packages/pressio_cmake_config.h.in pressio_cmake_config.h @ONLY)
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/pressio_cmake_config.h
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)

#=====================================================================
# Add uninstall target
#=====================================================================
add_custom_target (uninstall COMMAND ${CMAKE_COMMAND} -P
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/pressioUninstall.cmake)

#=====================================================================
# log to screen
#=====================================================================
message(STATUS "")
message(STATUS "${Green}Configure completed ${ColourReset}")
string(TOLOWER "${CMAKE_GENERATOR}" cmake_generator_tolower)
# if(cmake_generator_tolower MATCHES "makefile")
#   message(STATUS "Things you can do now:")
#   message(STATUS "----------------|--------------------------------------------------------")
#   message(STATUS "Command         |   Description")
#   message(STATUS "----------------|--------------------------------------------------------")
#   message(STATUS "make install    | Install headers to:")
#   message(STATUS "                |    ${CMAKE_INSTALL_PREFIX}/include")
#   message(STATUS "                | Change the install location using:")
#   message(STATUS "                |    cmake . -DCMAKE_INSTALL_PREFIX=yourprefix")
#   message(STATUS "                | ")
#   message(STATUS "make uninstall  | Removes files installed by make install")
#   message(STATUS "----------------|--------------------------------------------------------")
# endif()
# message(STATUS "")
