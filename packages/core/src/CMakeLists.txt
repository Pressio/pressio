
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}) #To include "core_config.h"

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------

###############
# HEADERS
###############
FILE(GLOB core_general_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
FILE(GLOB core_io_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/io/*.hpp")
FILE(GLOB core_shared_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/shared_base/*.hpp")

FILE(GLOB core_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/meta/*.hpp")

FILE(GLOB core_containers_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/*.hpp")
FILE(GLOB core_containers_dot_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/dot/*.hpp")
FILE(GLOB core_containers_dot_mvec_mvec_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/dot/mvec_mvec/*.hpp")
FILE(GLOB core_containers_dot_mvec_vec_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/dot/mvec_vec/*.hpp")
FILE(GLOB core_containers_dot_vec_mvec_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/dot/vec_mvec/*.hpp")
FILE(GLOB core_containers_dot_vec_vec_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/dot/vec_vec/*.hpp")
FILE(GLOB core_containers_dot_self_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/dot_self/*.hpp")
FILE(GLOB core_containers_dot_self_mvec_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/dot_self/mvec/*.hpp")
FILE(GLOB core_containers_mat_mat_prod_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/mat_mat_product/*.hpp")
FILE(GLOB core_containers_mat_vec_prod_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/mat_vec_product/*.hpp")
FILE(GLOB core_containers_multivec_vec_prod_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/multivec_vec_product/*.hpp")
FILE(GLOB core_containers_mat_multivec_prod_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/mat_multivec_product/*.hpp")
FILE(GLOB core_containers_norms_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/norms/*.hpp")
FILE(GLOB core_containers_minmax_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/min_max/*.hpp")

FILE(GLOB core_vector_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/vector/*.hpp")
FILE(GLOB core_vector_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/vector/meta/*.hpp")
FILE(GLOB core_vector_SHMEM_EXPR_T_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/vector/sharedmem_expr_templates/*.hpp")
FILE(GLOB core_vector_DISTMEM_EXPR_T_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/vector/distributed_expr_templates/*.hpp")
FILE(GLOB core_vector_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/vector/base/*.hpp")
FILE(GLOB core_vector_concrete_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/vector/concrete/*.hpp")

FILE(GLOB core_multi_vector_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/multi_vector/*.hpp")
FILE(GLOB core_multi_vector_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/multi_vector/meta/*.hpp")
FILE(GLOB core_multi_vector_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/multi_vector/base/*.hpp")
FILE(GLOB core_multi_vector_concrete_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/multi_vector/concrete/*.hpp")

FILE(GLOB core_matrix_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/matrix/*.hpp")
FILE(GLOB core_matrix_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/matrix/meta/*.hpp")
FILE(GLOB core_matrix_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/matrix/base/*.hpp")
FILE(GLOB core_matrix_concrete_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/matrix/concrete/*.hpp")

#################
#### SOURCES
#################
FILE(GLOB core_general_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc")
SET(SOURCES ${core_general_SRCS})

#--------------------------------------
SET(LUMPED_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/../../CORE_ALL
  ${CMAKE_CURRENT_SOURCE_DIR}/../../CORE_BASIC
  ${CMAKE_CURRENT_SOURCE_DIR}/../../CORE_MATRIX
  ${CMAKE_CURRENT_SOURCE_DIR}/../../CORE_MULTI_VECTOR
  ${CMAKE_CURRENT_SOURCE_DIR}/../../CORE_OPS
  ${CMAKE_CURRENT_SOURCE_DIR}/../../CORE_VECTOR
)
INSTALL(FILES ${LUMPED_HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

#--------------------------------------
# install the actual core headers
#--------------------------------------
SET(CORE_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/core/src)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h DESTINATION ${CORE_INCDIR})
INSTALL(FILES ${core_general_HEADERS} DESTINATION ${CORE_INCDIR})
INSTALL(FILES ${core_io_HEADERS} DESTINATION ${CORE_INCDIR}/io)
INSTALL(FILES ${core_shared_base_HEADERS} DESTINATION ${CORE_INCDIR}/shared_base)
INSTALL(FILES ${core_meta_HEADERS} DESTINATION ${CORE_INCDIR}/meta)

INSTALL(FILES ${core_containers_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops)
INSTALL(FILES ${core_containers_dot_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/dot)
INSTALL(FILES ${core_containers_dot_mvec_mvec_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/dot/mvec_mvec)
INSTALL(FILES ${core_containers_dot_mvec_vec_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/dot/mvec_vec)
INSTALL(FILES ${core_containers_dot_vec_mvec_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/dot/vec_mvec)
INSTALL(FILES ${core_containers_dot_vec_vec_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/dot/vec_vec)
INSTALL(FILES ${core_containers_dot_self_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/dot_self)
INSTALL(FILES ${core_containers_dot_self_mvec_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/dot_self/mvec)
INSTALL(FILES ${core_containers_mat_mat_prod_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/mat_mat_product)
INSTALL(FILES ${core_containers_mat_vec_prod_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/mat_vec_product)
INSTALL(FILES ${core_containers_multivec_vec_prod_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/multivec_vec_product)
INSTALL(FILES ${core_containers_mat_multivec_prod_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/mat_multivec_product)
INSTALL(FILES ${core_containers_norms_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/norms)
INSTALL(FILES ${core_containers_minmax_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/min_max)

INSTALL(FILES ${core_vector_HEADERS} DESTINATION ${CORE_INCDIR}/vector)
INSTALL(FILES ${core_vector_meta_HEADERS} DESTINATION ${CORE_INCDIR}/vector/meta)
INSTALL(FILES ${core_vector_SHMEM_EXPR_T_HEADERS} DESTINATION ${CORE_INCDIR}/vector/sharedmem_expr_templates)
INSTALL(FILES ${core_vector_DISTMEM_EXPR_T_HEADERS} DESTINATION ${CORE_INCDIR}/vector/distributed_expr_templates)
INSTALL(FILES ${core_vector_base_HEADERS} DESTINATION ${CORE_INCDIR}/vector/base)
INSTALL(FILES ${core_vector_concrete_HEADERS} DESTINATION ${CORE_INCDIR}/vector/concrete)

INSTALL(FILES ${core_multi_vector_HEADERS} DESTINATION ${CORE_INCDIR}/multi_vector)
INSTALL(FILES ${core_multi_vector_meta_HEADERS} DESTINATION ${CORE_INCDIR}/multi_vector/meta)
INSTALL(FILES ${core_multi_vector_base_HEADERS} DESTINATION ${CORE_INCDIR}/multi_vector/base)
INSTALL(FILES ${core_multi_vector_concrete_HEADERS} DESTINATION ${CORE_INCDIR}/multi_vector/concrete)

INSTALL(FILES ${core_matrix_HEADERS} DESTINATION ${CORE_INCDIR}/matrix)
INSTALL(FILES ${core_matrxi_meta_HEADERS} DESTINATION ${CORE_INCDIR}/matrix/meta)
INSTALL(FILES ${core_matrix_base_HEADERS} DESTINATION ${CORE_INCDIR}/matrix/base)
INSTALL(FILES ${core_matrix_concrete_HEADERS} DESTINATION ${CORE_INCDIR}/matrix/concrete)


#----------------------------------------------------------------
# C) Define the targets for package's library(s)
#----------------------------------------------------------------
TRIBITS_ADD_LIBRARY(
  core
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
