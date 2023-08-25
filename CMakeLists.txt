cmake_minimum_required(VERSION 3.18.0 FATAL_ERROR)
project(pressio C CXX)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include(colors)

#=====================================================================
# versioning
#=====================================================================
# adapted from Eigen
file(READ "${PROJECT_SOURCE_DIR}/include/pressio/macros.hpp" _pressio_macros)
string(REGEX MATCH "define[ \t]+PRESSIO_MAJOR_VERSION[ \t]+([0-9]+)" _pressio_major_version_match "${_pressio_macros}")
set(PRESSIO_MAJOR_VERSION "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+PRESSIO_MINOR_VERSION[ \t]+([0-9]+)" _pressio_minor_version_match "${_pressio_macros}")
set(PRESSIO_MINOR_VERSION "${CMAKE_MATCH_1}")
string(REGEX MATCH "define[ \t]+PRESSIO_PATCH_VERSION[ \t]+([0-9]+)" _pressio_patch_version_match "${_pressio_macros}")
set(PRESSIO_PATCH_VERSION "${CMAKE_MATCH_1}")
set(PRESSIO_VERSION_NUMBER ${PRESSIO_MAJOR_VERSION}.${PRESSIO_MINOR_VERSION}.${PRESSIO_PATCH_VERSION})
message("${Magenta}>> PRESSIO: version = ${PRESSIO_VERSION_NUMBER} ${ColourReset}")

#=====================================================================
# c++ standard
#=====================================================================
IF (NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
ENDIF()

message("${Magenta}>> PRESSIO: CMAKE_CXX_STANDARD is set to ${CMAKE_CXX_STANDARD} ${ColourReset}")

IF(CMAKE_CXX_STANDARD STREQUAL "14")
  if(NOT "cxx_std_14" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    message(FATAL_ERROR "Compiler does not support C++14. This is required.")
  endif()

ELSEIF(CMAKE_CXX_STANDARD STREQUAL "17")
  if(NOT "cxx_std_17" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    message(FATAL_ERROR "Compiler does not support C++17. This is required.")
  endif()

ELSEIF(CMAKE_CXX_STANDARD STREQUAL "20")
  message("${Magenta}>> PRESSIO: CMAKE_CXX_STANDARD = 20 ${ColourReset}")
  if(NOT "cxx_std_20" IN_LIST CMAKE_CXX_COMPILE_FEATURES)
    message(FATAL_ERROR "Compiler does not support C++20. This is required.")
  endif()
ENDIF()

################################################################################

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(
  DIRECTORY include/pressio
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(pressio INTERFACE)

target_include_directories(
  pressio INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(
  TARGETS pressio
  EXPORT pressioTargets
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(
  EXPORT pressioTargets
  FILE pressioTargets.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pressio)
export(
  TARGETS pressio
  FILE pressioTargets.cmake)

configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/pressioConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/pressioConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pressio)

write_basic_package_version_file (pressioConfigVersion.cmake
                                  VERSION ${PRESSIO_VERSION_NUMBER}
                                  COMPATIBILITY SameMajorVersion
				  ARCH_INDEPENDENT)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/pressioConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/pressioConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/pressio)

################################################################################

if(PRESSIO_ENABLE_FUNCTIONAL_TESTS OR PRESSIO_ENABLE_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif()
