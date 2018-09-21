
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------
#To include "core_config.h"
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/core)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

###############
# HEADERS
###############
FILE(GLOB core_general_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/*.hpp")
FILE(GLOB core_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/meta/*.hpp")
FILE(GLOB core_shared_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/shared_base/*.hpp")

FILE(GLOB core_containers_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/containers_ops/*.hpp")

FILE(GLOB core_vector_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/vector/*.hpp")
FILE(GLOB core_vector_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/vector/base/*.hpp")
FILE(GLOB core_vector_concrete_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/vector/concrete/*.hpp")

FILE(GLOB core_multi_vector_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/multi_vector/*.hpp")
FILE(GLOB core_multi_vector_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/multi_vector/base/*.hpp")
FILE(GLOB core_multi_vector_concrete_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/multi_vector/concrete/*.hpp")

FILE(GLOB core_matrix_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/matrix/*.hpp")
FILE(GLOB core_matrix_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/matrix/base/*.hpp")
FILE(GLOB core_matrix_concrete_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/matrix/concrete/*.hpp")

  
# SET(HEADERS
#   ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
#   ${core_general_HEADERS}
#   ${core_meta_HEADERS}
#   ${core_core_shared_base_HEADERS}
#   #
#   ${core_containers_ops_HEADERS}
#   #
#   ${core_vector_HEADERS}
#   ${core_vector_base_HEADERS}
#   ${core_vector_concrete_HEADERS}
#   #
#   ${core_multi_vector_HEADERS}
#   ${core_multi_vector_base_HEADERS}
#   ${core_multi_vector_concrete_HEADERS}
#   #
#   ${core_matrix_HEADERS}
#   ${core_matrix_base_HEADERS}
#   ${core_matrix_concrete_HEADERS}
# )


#################
#### SOURCES
#################
FILE(GLOB core_general_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/src/*.cc")
FILE(GLOB core_meta_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/meta/src/*.cc")
FILE(GLOB core_shared_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/shared_base/src/*.cc")

FILE(GLOB core_containers_ops_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/containers_ops_ops/src/*.cc")

FILE(GLOB core_vector_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/vector/src/*.cc")
FILE(GLOB core_vector_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/vector/base/src/*.cc")
FILE(GLOB core_vector_concrete_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/vector/concrete/src/*.cc")

FILE(GLOB core_multi_vector_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/multi_vector/src/*.cc")
FILE(GLOB core_multi_vector_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/multi_vector/base/src/*.cc")
FILE(GLOB core_multi_vector_concrete_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/multi_vector/concrete/src/*.cc")

FILE(GLOB core_matrix_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/matrix/src/*.cc")
FILE(GLOB core_matrix_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/matrix/base/src/*.cc")
FILE(GLOB core_matrix_concrete_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/core/matrix/concrete/src/*.cc")


SET(SOURCES 
  ${core_general_SRCS}
  ${core_meta_SRCS}
  ${core_core_shared_base_SRCS}
  #
  ${core_containers_ops_SRCS}
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


SET(CORE_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

#--------------------------------------
SET(LUMPED_HEADERS
  CORE_ALL
  CORE_BASIC
  CORE_MATRIX
  CORE_MULTI_VECTOR
  CORE_OPS
  CORE_VECTOR
)
INSTALL(FILES ${LUMPED_HEADERS} DESTINATION ${CORE_INCDIR})


#--------------------------------------
# install the actual core headers
#--------------------------------------
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h DESTINATION ${CORE_INCDIR}/core)

INSTALL(FILES ${core_general_HEADERS} DESTINATION ${CORE_INCDIR}/core)
INSTALL(FILES ${core_meta_HEADERS} DESTINATION ${CORE_INCDIR}/core/meta)
INSTALL(FILES ${core_core_shared_base_HEADERS} DESTINATION ${CORE_INCDIR}/core/shared_base)
  #
INSTALL(FILES ${core_containers_ops_HEADERS} DESTINATION ${CORE_INCDIR}/core/containers_ops)
  #
INSTALL(FILES ${core_vector_HEADERS} DESTINATION ${CORE_INCDIR}/core/vector)
INSTALL(FILES ${core_vector_base_HEADERS} DESTINATION ${CORE_INCDIR}/core/vector/base)
INSTALL(FILES ${core_vector_concrete_HEADERS} DESTINATION ${CORE_INCDIR}/core/vector/concrete)
  #
INSTALL(FILES ${core_multi_vector_HEADERS} DESTINATION ${CORE_INCDIR}/core/multi_vector)
INSTALL(FILES ${core_multi_vector_base_HEADERS} DESTINATION ${CORE_INCDIR}/core/multi_vector/base)
INSTALL(FILES ${core_multi_vector_concrete_HEADERS} DESTINATION ${CORE_INCDIR}/core/multi_vector/concrete)
  #
INSTALL(FILES ${core_matrix_HEADERS} DESTINATION ${CORE_INCDIR}/core/matrix)
INSTALL(FILES ${core_matrix_base_HEADERS} DESTINATION ${CORE_INCDIR}/core/matrix/base)
INSTALL(FILES ${core_matrix_concrete_HEADERS} DESTINATION ${CORE_INCDIR}/core/matrix/concrete)


TRIBITS_ADD_LIBRARY(
  core
  SOURCES ${SOURCES}
  DEPLIBS
  )


# #----------------------------------------------------------------
# # C) Define the targets for package's library(s)
# #----------------------------------------------------------------
# TRIBITS_ADD_LIBRARY(
#   core
#   HEADERS ${HEADERS}
#   SOURCES ${SOURCES}
#   )

# #
# # D) Define an executable that gets installed
# #
# TRIBITS_ADD_EXECUTABLE(commomBase-helloworld  NOEXEPREFIX  NOEXESUFFIX
#   SOURCES  SimpleCxx_HelloWorld_Prog.cpp
#   INSTALLABLE
#   )
