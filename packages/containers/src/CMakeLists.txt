
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}) #To include "containers_config.h"

# Set the base directory in order to find pressio_version.h
#ASSERT_DEFINED(PROJECT_SOURCE_DIR)
#INCLUDE_DIRECTORIES(${${PROJECT_NAME}_BINARY_DIR})

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------

###############
# HEADERS
###############
FILE(GLOB containers_general_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
FILE(GLOB containers_shared_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/shared_base/*.hpp")

FILE(GLOB containers_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/meta/*.hpp")

FILE(GLOB containers_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/ops/*.hpp")

FILE(GLOB containers_eigen_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/ops/eigen/*.hpp")
FILE(GLOB containers_epetra_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/ops/epetra/*.hpp")
FILE(GLOB containers_tpetra_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/ops/tpetra/*.hpp")
FILE(GLOB containers_tpetra_block_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/ops/tpetra_block/*.hpp")
FILE(GLOB containers_teuchos_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/ops/teuchos/*.hpp")
FILE(GLOB containers_kokkos_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/ops/kokkos/*.hpp")
FILE(GLOB containers_pybind11_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/ops/pybind11/*.hpp")

FILE(GLOB containers_vector_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/vector/*.hpp")
FILE(GLOB containers_vector_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/vector/meta/*.hpp")
FILE(GLOB containers_vector_SHMEM_EXPR_T_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/vector/sharedmem_expr_templates/*.hpp")
FILE(GLOB containers_vector_DISTMEM_EXPR_T_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/vector/distributed_expr_templates/*.hpp")
FILE(GLOB containers_vector_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/vector/base/*.hpp")
FILE(GLOB containers_vector_concrete_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/vector/concrete/*.hpp")

FILE(GLOB containers_multi_vector_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/multi_vector/*.hpp")
FILE(GLOB containers_multi_vector_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/multi_vector/meta/*.hpp")
FILE(GLOB containers_multi_vector_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/multi_vector/base/*.hpp")
FILE(GLOB containers_multi_vector_concrete_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/multi_vector/concrete/*.hpp")

FILE(GLOB containers_matrix_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/matrix/*.hpp")
FILE(GLOB containers_matrix_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/matrix/meta/*.hpp")
FILE(GLOB containers_matrix_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/matrix/base/*.hpp")
FILE(GLOB containers_matrix_concrete_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/matrix/concrete/*.hpp")

#################
#### SOURCES
#################
FILE(GLOB containers_general_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")
SET(SOURCES ${containers_general_SRCS})

#--------------------------------------
SET(LUMPED_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/../../CONTAINERS_ALL
  ${CMAKE_CURRENT_SOURCE_DIR}/../../CONTAINERS_BASIC
  ${CMAKE_CURRENT_SOURCE_DIR}/../../CONTAINERS_MATRIX
  ${CMAKE_CURRENT_SOURCE_DIR}/../../CONTAINERS_MULTI_VECTOR
  ${CMAKE_CURRENT_SOURCE_DIR}/../../CONTAINERS_OPS
  ${CMAKE_CURRENT_SOURCE_DIR}/../../CONTAINERS_VECTOR
)
INSTALL(FILES ${LUMPED_HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

#--------------------------------------
# install the actual containers headers
#--------------------------------------
SET(CONTAINERS_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/containers/src)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h DESTINATION ${CONTAINERS_INCDIR})
INSTALL(FILES ${containers_general_HEADERS} DESTINATION ${CONTAINERS_INCDIR})
INSTALL(FILES ${containers_shared_base_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/shared_base)
INSTALL(FILES ${containers_meta_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/meta)

INSTALL(FILES ${containers_ops_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/ops)
INSTALL(FILES ${containers_eigen_ops_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/ops/eigen)
INSTALL(FILES ${containers_epetra_ops_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/ops/epetra)
INSTALL(FILES ${containers_tpetra_ops_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/ops/tpetra)
INSTALL(FILES ${containers_tpetra_block_ops_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/ops/tpetra_block)
INSTALL(FILES ${containers_teuchos_ops_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/ops/teuchos)
INSTALL(FILES ${containers_kokkos_ops_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/ops/kokkos)
INSTALL(FILES ${containers_pybind11_ops_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/ops/pybind11)

INSTALL(FILES ${containers_vector_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/vector)
INSTALL(FILES ${containers_vector_meta_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/vector/meta)
INSTALL(FILES ${containers_vector_SHMEM_EXPR_T_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/vector/sharedmem_expr_templates)
INSTALL(FILES ${containers_vector_DISTMEM_EXPR_T_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/vector/distributed_expr_templates)
INSTALL(FILES ${containers_vector_base_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/vector/base)
INSTALL(FILES ${containers_vector_concrete_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/vector/concrete)

INSTALL(FILES ${containers_multi_vector_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/multi_vector)
INSTALL(FILES ${containers_multi_vector_meta_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/multi_vector/meta)
INSTALL(FILES ${containers_multi_vector_base_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/multi_vector/base)
INSTALL(FILES ${containers_multi_vector_concrete_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/multi_vector/concrete)

INSTALL(FILES ${containers_matrix_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/matrix)
INSTALL(FILES ${containers_matrix_meta_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/matrix/meta)
INSTALL(FILES ${containers_matrix_base_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/matrix/base)
INSTALL(FILES ${containers_matrix_concrete_HEADERS} DESTINATION ${CONTAINERS_INCDIR}/matrix/concrete)

#----------------------------------------------------------------
# C) Define the targets for package's library(s)
#----------------------------------------------------------------
TRIBITS_ADD_LIBRARY(
  containers
  #   HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS
  )

# #
# # D) Define an executable that gets installed
# #
# TRIBITS_ADD_EXECUTABLE(commomBase-helloworld  NOEXEPREFIX  NOEXESUFFIX
#   SOURCES  SimpleCxx_HelloWorld_Prog.cpp
#   INSTALLABLE
#   )
