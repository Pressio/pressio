
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
FILE(GLOB core_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/meta/*.hpp")
FILE(GLOB core_tinympl_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/meta/tinympl/*.hpp")
FILE(GLOB core_tinymplvar_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/meta/tinympl/variadic/*.hpp")

FILE(GLOB core_shared_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/shared_base/*.hpp")

FILE(GLOB core_containers_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/*.hpp")
FILE(GLOB core_containers_dot_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/dot/*.hpp")
FILE(GLOB core_containers_dot_mvec_mvec_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/dot/mvec_mvec/*.hpp")
FILE(GLOB core_containers_dot_mvec_vec_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/dot/mvec_vec/*.hpp")

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
FILE(GLOB core_vector_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/vector/base/*.hpp")
FILE(GLOB core_vector_concrete_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/vector/concrete/*.hpp")

FILE(GLOB core_multi_vector_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/multi_vector/*.hpp")
FILE(GLOB core_multi_vector_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/multi_vector/base/*.hpp")
FILE(GLOB core_multi_vector_concrete_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/multi_vector/concrete/*.hpp")

FILE(GLOB core_matrix_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/matrix/*.hpp")
FILE(GLOB core_matrix_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/matrix/base/*.hpp")
FILE(GLOB core_matrix_concrete_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/matrix/concrete/*.hpp")


#################
#### SOURCES
#################
FILE(GLOB core_general_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc")
FILE(GLOB core_meta_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/meta/src/*.cc")
FILE(GLOB core_shared_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/shared_base/src/*.cc")

FILE(GLOB core_containers_ops_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/src/*.cc")
FILE(GLOB core_containers_dot_ops_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/dot/src/*.cc")
FILE(GLOB core_containers_dot_mvec_mvec_ops_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/dot/mvec_mvec//src/*.cc")
FILE(GLOB core_containers_dot_mvec_vec_ops_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/dot/mvec_vec/src/*.cc")

FILE(GLOB core_containers_mat_mat_prod_ops_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/mat_mat_product/src/*.cc")
FILE(GLOB core_containers_mat_vec_prod_ops_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/mat_vec_product/src/*.cc")
FILE(GLOB core_containers_multivec_vec_prod_ops_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/multivec_vec_product/src/*.cc")
FILE(GLOB core_containers_mat_multivec_prod_ops_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/mat_multivec_product/src/*.cc")
FILE(GLOB core_containers_norms_ops_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/norms/src/*.cc")
FILE(GLOB core_containers_minmax_ops_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/containers_ops/min_max/src/*.cc")


FILE(GLOB core_vector_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/vector/src/*.cc")
FILE(GLOB core_vector_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/vector/base/src/*.cc")
FILE(GLOB core_vector_concrete_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/vector/concrete/src/*.cc")

FILE(GLOB core_multi_vector_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/multi_vector/src/*.cc")
FILE(GLOB core_multi_vector_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/multi_vector/base/src/*.cc")
FILE(GLOB core_multi_vector_concrete_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/multi_vector/concrete/src/*.cc")

FILE(GLOB core_matrix_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/matrix/src/*.cc")
FILE(GLOB core_matrix_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/matrix/base/src/*.cc")
FILE(GLOB core_matrix_concrete_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/matrix/concrete/src/*.cc")


SET(SOURCES
  ${core_general_SRCS}
  ${core_meta_SRCS}
  ${core_core_shared_base_SRCS}
  #
  ${core_containers_ops_SRCS}
  ${core_containers_dot_ops_SRCS}
  ${core_containers_dot_mvec_mvec_ops_SRCS}
  ${core_containers_dot_mvec_vec_ops_SRCS}
  ${core_containers_mat_mat_prod_ops_SRCS}
  ${core_containers_mat_vec_prod_ops_SRCS}
  ${core_containers_multivec_vec_prod_ops_SRCS}
  ${core_containers_mat_multivec_prod_ops_SRCS}
  ${core_containers_norms_ops_SRCS}
  ${core_containers_minmax_ops_SRCS}
  #
  ${core_vector_SRCS}
  ${core_vector_base_SRCS}
  ${core_vector_concrete_SRCS}
  #
  ${core_multi_vector_SRCS}
  ${core_multi_vector_base_SRCS}
  ${core_multi_vector_concrete_SRCS}
  #
  ${core_matrix_SRCS}
  ${core_matrix_base_SRCS}
  ${core_matrix_concrete_SRCS}
)


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
INSTALL(FILES ${core_meta_HEADERS} DESTINATION ${CORE_INCDIR}/meta)
INSTALL(FILES ${core_tinympl_HEADERS} DESTINATION ${CORE_INCDIR}/meta/tinympl)
INSTALL(FILES ${core_tinymplvar_HEADERS} DESTINATION ${CORE_INCDIR}/meta/tinympl/variadic)

INSTALL(FILES ${core_shared_base_HEADERS} DESTINATION ${CORE_INCDIR}/shared_base)
#
INSTALL(FILES ${core_containers_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops)

INSTALL(FILES ${core_containers_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops)
INSTALL(FILES ${core_containers_dot_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/dot)
INSTALL(FILES ${core_containers_dot_mvec_mvec_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/dot/mvec_mvec)
INSTALL(FILES ${core_containers_dot_mvec_vec_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/dot/mvec_vec)

INSTALL(FILES ${core_containers_mat_mat_prod_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/mat_mat_product)
INSTALL(FILES ${core_containers_mat_vec_prod_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/mat_vec_product)
INSTALL(FILES ${core_containers_multivec_vec_prod_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/multivec_vec_product)
INSTALL(FILES ${core_containers_mat_multivec_prod_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/mat_multivec_product)
INSTALL(FILES ${core_containers_norms_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/norms)
INSTALL(FILES ${core_containers_minmax_ops_HEADERS} DESTINATION ${CORE_INCDIR}/containers_ops/min_max)

#
INSTALL(FILES ${core_vector_HEADERS} DESTINATION ${CORE_INCDIR}/vector)
INSTALL(FILES ${core_vector_base_HEADERS} DESTINATION ${CORE_INCDIR}/vector/base)
INSTALL(FILES ${core_vector_concrete_HEADERS} DESTINATION ${CORE_INCDIR}/vector/concrete)
#
INSTALL(FILES ${core_multi_vector_HEADERS} DESTINATION ${CORE_INCDIR}/multi_vector)
INSTALL(FILES ${core_multi_vector_base_HEADERS} DESTINATION ${CORE_INCDIR}/multi_vector/base)
INSTALL(FILES ${core_multi_vector_concrete_HEADERS} DESTINATION ${CORE_INCDIR}/multi_vector/concrete)
#
INSTALL(FILES ${core_matrix_HEADERS} DESTINATION ${CORE_INCDIR}/matrix)
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
