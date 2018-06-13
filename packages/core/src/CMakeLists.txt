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
  core_forward_declarations.hpp
  core_static_assert_definitions.hpp
  core_crtp_helper.hpp
  core_operators_base.hpp
  meta/core_meta_basic.hpp
  meta/core_meta_detect_operators.hpp
  meta/core_meta_detect_typedefs.hpp
  )

APPEND_SET(HEADERS
  vector/core_vector_traits.hpp
  vector/core_vector_meta.hpp  
  vector/base/core_vector_generic_base.hpp
  vector/base/core_vector_math_base.hpp
  vector/base/core_vector_serial_base.hpp
  vector/core_vector_distributed_epetra.hpp
  vector/core_vector_serial_eigen.hpp
  vector/core_vector_serial_stdlib.hpp
  vector/core_vector_serial_userdefined.hpp
  vector/core_vector_dot_product.hpp
  )

APPEND_SET(HEADERS
  matrix/core_matrix_traits.hpp
  matrix/core_matrix_meta.hpp  
  matrix/base/core_matrix_dense_distributed_base.hpp
  matrix/base/core_matrix_dense_serial_base.hpp
  matrix/base/core_matrix_generic_base.hpp
  matrix/base/core_matrix_sparse_distributed_base.hpp
  matrix/base/core_matrix_sparse_serial_base.hpp
  matrix/core_matrix_dense_serial_eigen.hpp
  matrix/core_matrix_dense_serial_stdlib.hpp
  matrix/core_matrix_sparse_serial_eigen.hpp
  matrix/core_matrix_vector_product.hpp
)



# adding all sources is temporary, we need to wait 
# until tribits allows one to have a header-only library.
APPEND_SET(SOURCES
  src/core_meta_basic.cc
  src/core_meta_detect_operators.cc
  src/core_meta_detect_typedefs.cc
  src/core_forward_declarations.cc
  src/core_static_assert_definitions.cc
  src/core_crtp_helper.cc
  src/core_operators_base.cc
  #
  vector/src/core_vector_traits.cc
  vector/src/core_vector_meta.cc
  vector/src/core_vector_generic_base.cc
  vector/src/core_vector_math_base.cc
  vector/src/core_vector_serial_base.cc
  vector/src/core_vector_distributed_epetra.cc
  vector/src/core_vector_serial_eigen.cc
  vector/src/core_vector_serial_stdlib.cc
  vector/src/core_vector_serial_userdefined.cc
  vector/src/core_vector_dot_product.cc
  #
  matrix/src/core_matrix_traits.cc
  matrix/src/core_matrix_meta.cc
  matrix/src/core_matrix_dense_distributed_base.cc
  matrix/src/core_matrix_dense_serial_base.cc
  matrix/src/core_matrix_generic_base.cc
  matrix/src/core_matrix_sparse_distributed_base.cc
  matrix/src/core_matrix_sparse_serial_base.cc
  matrix/src/core_matrix_dense_serial_eigen.cc
  matrix/src/core_matrix_dense_serial_stdlib.cc
  matrix/src/core_matrix_sparse_serial_eigen.cc
  matrix/src/core_matrix_vector_product.cc
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
