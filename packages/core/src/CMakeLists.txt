#
# A) Package-specific configuration options
#
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and include directories)
#
SET(HEADERS "")
SET(SOURCES "")

# To include the configured header above
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  core_ConfigDefs.hpp
  core_meta.hpp
  core_forward_declarations.hpp
  core_static_assert.hpp
  #
  vector/core_vector_traits.hpp
  vector/core_vector_generic_base.hpp
  vector/core_vector_dist_base.hpp
  vector/core_vector_serial_base.hpp
  vector/core_vector_math_base.hpp
  #
  vector/core_vector_epetra.hpp
  vector/core_vector_std.hpp
  vector/core_vector_serial_arbitrary.hpp
  vector/core_vector_eigen.hpp
  #
  matrix/core_matrix_traits.hpp
  matrix/core_matrix_generic_base.hpp
  matrix/core_matrix_eigen.hpp
)

APPEND_SET(SOURCES
  src/core_meta.cc
  src/core_forward_declarations.cc
  src/core_static_assert.cc
  #
  vector/src/core_vector_traits.cc
  vector/src/core_vector_generic_base.cc
  vector/src/core_vector_dist_base.cc
  vector/src/core_vector_serial_base.cc
  vector/src/core_vector_math_base.cc
  #
  vector/src/core_vector_epetra.cc
  vector/src/core_vector_std.cc
  vector/src/core_vector_serial_arbitrary.cc
  vector/src/core_vector_eigen.cc
  #
  matrix/src/core_matrix_traits.cc
  matrix/src/core_matrix_generic_base.cc
  matrix/src/core_matrix_eigen.cc
)

#
# C) Define the targets for package's library(s)
#

# # # this seems to work for header only library
# TRIBITS_INSTALL_HEADERS(
#   HEADERS ${HEADERS}
# )
# this is the original in case we have sources too for the library
TRIBITS_ADD_LIBRARY(
  core
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

# #
# # D) Define an executable that gets installed
# #
# TRIBITS_ADD_EXECUTABLE(commomBase-helloworld  NOEXEPREFIX  NOEXESUFFIX
#   SOURCES  SimpleCxx_HelloWorld_Prog.cpp
#   INSTALLABLE
#   )
