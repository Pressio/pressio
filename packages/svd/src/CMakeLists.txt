
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}) #To include "svd_config.h"

# Set the base directory in order to find pressio_version.h
ASSERT_DEFINED(PROJECT_SOURCE_DIR)
INCLUDE_DIRECTORIES(${PROJECT_NAME}_BINARY_DIR})

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------

###############
# HEADERS
###############
FILE(GLOB svd_general_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

#################
#### SOURCES
#################
FILE(GLOB svd_general_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")
SET(SOURCES 
  ${svd_general_SRCS}
)


#--------------------------------------
SET(LUMPED_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/../../SVD_BASIC
)
INSTALL(FILES ${LUMPED_HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

#--------------------------------------
# install the actual headers
#--------------------------------------
SET(SVD_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/svd/src)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h DESTINATION ${SVD_INCDIR})
INSTALL(FILES ${svd_general_HEADERS} DESTINATION ${SVD_INCDIR})

#----------------------------------------------------------------
# C) Define the targets for package's library(s)
#----------------------------------------------------------------
TRIBITS_ADD_LIBRARY(
  svd
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
