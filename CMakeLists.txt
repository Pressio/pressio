
cmake_minimum_required(VERSION 3.11.0 FATAL_ERROR)
project(pressio C CXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake;${PROJECT_SOURCE_DIR}/cmake/tpls)

#=====================================================================
# do some guarding here for specific cases
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
if(    NOT cmake_build_type_tolower STREQUAL "debug"
   AND NOT cmake_build_type_tolower STREQUAL "release")
  message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\".
    Allowed values are Debug, Release (case-insensitive).")
endif()


#=====================================================================
# we need c++11
#=====================================================================
include(CheckCXXCompilerFlag)

check_cxx_compiler_flag("-std=c++11" COMPILER_SUPPORT_CPP11)
if(NOT COMPILER_SUPPORT_CPP11)
  message(FATAL_ERROR "Compiler does not support -std=c++11. This is required.")
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


#=====================================================================
# other stuff
#=====================================================================
include(GNUInstallDirs)
include_directories(${CMAKE_CURRENT_BINARY_DIR})
add_subdirectory(packages)


#=====================================================================
# some optional and TPLs
# note that linking to TPLs is only needed when building tests
#=====================================================================
if( cmake_build_type_tolower STREQUAL "debug" AND NOT DEBUG_PRINT)
  message( "CMAKE_BUILD_TYPE==Debug, setting DEBUG_PRINT=ON" )
  message( "If you want to turn off debug prints, reconfigure with -DDEBUG_PRINT=OFF" )
  set(DEBUG_PRINT ON)
endif()

option(BUILD_UNIT_TESTS "Enabling pressio unit testing." OFF)
option(BUILD_TESTS "Enable ALL pressio tests: regression and unit tests." OFF)

# ALL_INC/LIB_DIR contain all include/lib directories as we add TPLs
# we fill these vars up based on which TPLs are enabled. 
# ALL_INC/LIB_DIR are then used for creating target executables
set(ALL_INC_DIR "")
set(ALL_LINK_DIR "")
set(ALL_LINK_LIB "")
include(tplEigen)
include(tplMPI)
include(tplTrilinos)
include(tplBlas)
include(tplLapack)
include(tplKokkos)
include(tplPybind11)
include(tplArmadillo)
include(tplBlaze)


#=====================================================================
# if we want to enable unit tests
#=====================================================================
if(BUILD_UNIT_TESTS OR BUILD_TESTS)
  # tell CMake to use CTest extension
  enable_testing()
  # add subdir where the unit tests are
  add_subdirectory(unit_tests)
endif()

#=====================================================================
# configure file for cmakedefines
#=====================================================================
configure_file(pressio_cmake_config.h.in pressio_cmake_config.h @ONLY)
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
message(STATUS "Configured pressio")
message(STATUS "")

string(TOLOWER "${CMAKE_GENERATOR}" cmake_generator_tolower)
if(cmake_generator_tolower MATCHES "makefile")
  message(STATUS "Things you can do now:")
  message(STATUS "----------------|--------------------------------------------------------")
  message(STATUS "Command         |   Description")
  message(STATUS "----------------|--------------------------------------------------------")
  message(STATUS "make install    | Install headers to:")
  message(STATUS "                |    ${CMAKE_INSTALL_PREFIX}/include")
  message(STATUS "                | Change the install location using:")
  message(STATUS "                |    cmake . -DCMAKE_INSTALL_PREFIX=yourprefix")
  message(STATUS "                | ")
  message(STATUS "make uninstall  | Removes files installed by make install")
  message(STATUS "----------------|--------------------------------------------------------")
endif()
message(STATUS "")
