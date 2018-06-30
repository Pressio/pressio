#
# A) Package-specific configuration options
#
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and include directories)
#
SET(HEADERS "")
SET(SOURCES "")

#To include "core_config.h"
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

###############
# HEADERS
###############
FILE(GLOB core_general_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
FILE(GLOB core_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/meta/*.hpp")
FILE(GLOB core_utils_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils/*.hpp")

FILE(GLOB core_vector_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/vector/*.hpp")
FILE(GLOB core_vector_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/vector/meta/*.hpp")
FILE(GLOB core_vector_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/vector/operations/*.hpp")
FILE(GLOB core_vector_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/vector/base/*.hpp")
FILE(GLOB core_vector_concrete_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/vector/concrete/*.hpp")

FILE(GLOB core_matrix_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/matrix/*.hpp")
FILE(GLOB core_matrix_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/matrix/meta/*.hpp")
FILE(GLOB core_matrix_ops_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/matrix/operations/*.hpp")
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
FILE(GLOB core_utils_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/utils/src/*.cc")

FILE(GLOB core_vector_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/vector/src/*.cc")
FILE(GLOB core_vector_meta_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/vector/meta/src/*.cc")
FILE(GLOB core_vector_ops_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/vector/operations/src/*.cc")
FILE(GLOB core_vector_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/vector/base/src/*.cc")
FILE(GLOB core_vector_concrete_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/vector/concrete/src/*.cc")


# MESSAGE( "GIGI" ${general_HEADERS} )
APPEND_SET(HEADERS
  ${core_general_HEADERS}
  ${core_meta_HEADERS}
  ${core_utils_HEADERS}
  ${core_vector_HEADERS}
  ${core_vector_meta_HEADERS}
  ${core_vector_ops_HEADERS}
  ${core_vector_base_HEADERS}
  ${core_vector_concrete_HEADERS}
  ${core_matrix_HEADERS}
  ${core_matrix_meta_HEADERS}
  ${core_matrix_ops_HEADERS}
  ${core_matrix_base_HEADERS}
  ${core_matrix_concrete_HEADERS}
)
APPEND_SET(SOURCES
  ${core_general_SRCS}
  ${core_meta_SRCS}
  ${core_utils_SRCS}
  ${core_vector_SRCS}
  ${core_vector_meta_SRCS}
  ${core_vector_ops_SRCS}
  ${core_vector_base_SRCS}
  ${core_vector_concrete_SRCS}
  ${core_matrix_SRCS}
  ${core_matrix_meta_SRCS}
  ${core_matrix_ops_SRCS}
  ${core_matrix_base_SRCS}
  ${core_matrix_concrete_SRCS}
)



# C) Define the targets for package's library(s)

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
