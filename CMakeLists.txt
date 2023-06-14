
cmake_minimum_required(VERSION 3.18.0 FATAL_ERROR)
project(pressio C CXX)

# setting where other cmake commands are
set(CMAKE_MODULE_PATH
  ${PROJECT_SOURCE_DIR}/cmake;${PROJECT_SOURCE_DIR}/cmake/tpls)

# include commands for coloring prints
include(colors)

#=====================================================================
# versioning
#=====================================================================
file(READ "version.txt" tmpver)
string(REGEX MATCH "VERSION_MAJOR ([0-9]*)" _ ${tmpver})
set(Pressio_VERSION_MAJOR ${CMAKE_MATCH_1})
string(REGEX MATCH "VERSION_MINOR ([0-9]*)" _ ${tmpver})
set(Pressio_VERSION_MINOR ${CMAKE_MATCH_1})
string(REGEX MATCH "VERSION_PATCH ([0-9]*)" _ ${tmpver})
set(Pressio_VERSION_PATCH ${CMAKE_MATCH_1})

set(Pressio_VERSION "${Pressio_VERSION_MAJOR}.${Pressio_VERSION_MINOR}.${Pressio_VERSION_PATCH}")
message("${Magenta}>> Pressio version = ${Pressio_VERSION} ${ColourReset}")

#=====================================================================
# guard against in-source builds
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
  message(FATAL_ERROR "In-source builds are not allowed.
    Please make a new directory (called a build directory) and run CMake from there.
    You may need to remove CMakeCache.txt.")
endif()

#=====================================================================
# check build type
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()

# convert cmake build type to lower string
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_tolower)
if( NOT cmake_build_type_tolower STREQUAL "debug"
    AND NOT cmake_build_type_tolower STREQUAL "release")
  message(FATAL_ERROR "Unknown build type \"${CMAKE_BUILD_TYPE}\".
    Allowed values are Debug, Release (case-insensitive).")
endif()

#=====================================================================
# c++ standard
#=====================================================================
SET(PRESSIO_ENABLE_CXX14 OFF)
SET(PRESSIO_ENABLE_CXX17 OFF)
SET(PRESSIO_ENABLE_CXX20 OFF)

IF (NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
ENDIF()

include(CheckCXXCompilerFlag)
IF(CMAKE_CXX_STANDARD STREQUAL "14")
  SET(PRESSIO_ENABLE_CXX14 ON)
  check_cxx_compiler_flag("-std=c++14" COMPILER_SUPPORT_TARGET_CPP14)
  if(NOT COMPILER_SUPPORT_TARGET_CPP14)
    message(FATAL_ERROR "Compiler does not support -std=c++14. This is required.")
  endif()
  MESSAGE(STATUS "Setting PRESSIO_ENABLE_CXX14=ON")

ELSEIF(CMAKE_CXX_STANDARD STREQUAL "17")
  SET(PRESSIO_ENABLE_CXX14 ON)
  SET(PRESSIO_ENABLE_CXX17 ON)
  check_cxx_compiler_flag("-std=c++17" COMPILER_SUPPORT_TARGET_CPP17)
  if(NOT COMPILER_SUPPORT_TARGET_CPP17)
    message(FATAL_ERROR "Compiler does not support -std=c++17. This is required.")
  endif()
  MESSAGE(STATUS "Setting PRESSIO_ENABLE_CXX14=ON")
  MESSAGE(STATUS "Setting PRESSIO_ENABLE_CXX17=ON")

ELSEIF(CMAKE_CXX_STANDARD STREQUAL "20")
  SET(PRESSIO_ENABLE_CXX14 ON)
  SET(PRESSIO_ENABLE_CXX17 ON)
  SET(PRESSIO_ENABLE_CXX20 ON)
  check_cxx_compiler_flag("-std=c++20" COMPILER_SUPPORT_TARGET_CPP20)
  if(NOT COMPILER_SUPPORT_TARGET_CPP20)
    message(FATAL_ERROR "Compiler does not support -std=c++20. This is required.")
  endif()
  MESSAGE(STATUS "Setting PRESSIO_ENABLE_CXX14=ON")
  MESSAGE(STATUS "Setting PRESSIO_ENABLE_CXX17=ON")
  MESSAGE(STATUS "Setting PRESSIO_ENABLE_CXX20=ON")

ENDIF()

MESSAGE("${Magenta}>> Setting CXX standard to ${CMAKE_CXX_STANDARD} ${ColourReset}")

#=====================================================================
# do processing
#=====================================================================
include(GNUInstallDirs)

message("")
message("${Magenta}>> Processing:${ColourReset}")

if( cmake_build_type_tolower STREQUAL "debug")
  set(PRESSIO_DEBUG_BUILD ON)
  message( ">> PRESSIO_DEBUG_BUILD is ON since CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}" )
  message("")
endif()

if( cmake_build_type_tolower STREQUAL "debug" AND NOT PRESSIO_ENABLE_DEBUG_PRINT)
  set(PRESSIO_ENABLE_DEBUG_PRINT ON)
  message( ">> CMAKE_BUILD_TYPE==Debug, setting PRESSIO_ENABLE_DEBUG_PRINT=${PRESSIO_ENABLE_DEBUG_PRINT}" )
  message( "  To disable: reconfigure with -DPRESSIO_ENABLE_DEBUG_PRINT=OFF" )
  message("")
endif()
if(PRESSIO_ENABLE_DEBUG_PRINT)
  message( ">> Enabling debug prints since PRESSIO_ENABLE_DEBUG_PRINT=ON" )
  message( "  To disable: reconfigure with -DPRESSIO_ENABLE_DEBUG_PRINT=OFF" )
  message("")
endif()

option(PRESSIO_ENABLE_TESTS "Enable ALL pressio tests." OFF)
option(PRESSIO_ENABLE_FUNCTIONAL_SMALL_TESTS "Enabling pressio functional small testing." OFF)
option(PRESSIO_ENABLE_FUNCTIONAL_MEDIUM_TESTS "Enabling pressio functional medium testing." OFF)
option(PRESSIO_ENABLE_FUNCTIONAL_LARGE_TESTS "Enabling pressio functional large testing." OFF)

#=====================================================================
# note that here we don't link anything, the TPLs are only needed
# when building tests or when the user builds their app.
# For just installing pressio, we only need to set the cmake vars.
#=====================================================================
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
# install
#=====================================================================
set(INC_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include/pressio)

configure_file(include/pressio/pressio_cmake_config.h.in pressio_cmake_config.h @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/pressio_cmake_config.h DESTINATION ${INC_INSTALL_DIR})

file(GLOB include_helpers "${CMAKE_CURRENT_SOURCE_DIR}/include/pressio/*.hpp")
install(FILES ${include_helpers} DESTINATION ${INC_INSTALL_DIR})

foreach(X IN ITEMS mpl utils type_traits concepts expressions ops qr solvers_linear solvers_nonlinear ode rom)
  install(
    DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/pressio/${X}
    DESTINATION ${INC_INSTALL_DIR} FILES_MATCHING PATTERN "*.hpp")
endforeach()

#=====================================================================
# Add uninstall target
#=====================================================================
add_custom_target (uninstall COMMAND ${CMAKE_COMMAND} -P
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/pressioUninstall.cmake)

if(PRESSIO_ENABLE_FUNCTIONAL_TESTS OR PRESSIO_ENABLE_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()

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
