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
  meta.hpp
  forwardDeclarations.hpp
  staticAssert.hpp
  #
  vector/core_vectorTraits.hpp
  vector/core_vectorBaseImpl.hpp
  vector/core_vectorDistImpl.hpp
  vector/core_vectorSerImpl.hpp
  vector/core_vectorMathImpl.hpp
  vector/core_vector_epetra.hpp
  vector/core_vector_std.hpp
  vector/core_vector_serarb.hpp
  vector/core_vector_eigen.hpp
  #
  # matrix/core_matrixTraits.hpp
  # matrix/core_matrixBaseImpl.hpp
  # matrix/core_matrix_eigen.hpp
  # matrix/core_matrix_eigenStatic.hpp
  # matrix/core_matrix_std.hpp
)

APPEND_SET(SOURCES
  src/meta.cc
  src/forwardDeclarations.cc
  src/staticAssert.cc
  #
  vector/src/core_vectorTraits.cc
  vector/src/core_vectorBaseImpl.cc
  vector/src/core_vectorDistImpl.cc
  vector/src/core_vectorSerImpl.cc
  vector/src/core_vectorMathImpl.cc
  vector/src/core_vector_epetra.cc
  vector/src/core_vector_std.cc
  vector/src/core_vector_serarb.cc
  vector/src/core_vector_eigen.cc
  #
  # matrix/src/core_matrixTraits.cc
  # matrix/src/core_matrixBaseImpl.cc
  # matrix/src/core_matrix_eigen.cc
  # matrix/src/core_matrix_eigenStatic.cc
  # matrix/src/core_matrix_std.cc
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
