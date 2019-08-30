
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}) #To include "utils_config.h"

## Set the base directory in order to find pressio_version.h
#ASSERT_DEFINED(PROJECT_SOURCE_DIR)
#INCLUDE_DIRECTORIES(${${PROJECT_NAME}_BINARY_DIR})

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------

###############
# HEADERS
###############
FILE(GLOB general_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
FILE(GLOB io_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/io/*.hpp")

#################
#### SOURCES
#################
FILE(GLOB general_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")
SET(SOURCES ${general_SRCS})

#--------------------------------------
SET(LUMPED_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/../../UTILS_ALL
)
INSTALL(FILES ${LUMPED_HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

#--------------------------------------
# install the actual headers
#--------------------------------------
SET(UTILS_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/utils/src)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h DESTINATION ${UTILS_INCDIR})
INSTALL(FILES ${general_HEADERS} DESTINATION ${UTILS_INCDIR})
INSTALL(FILES ${io_HEADERS} DESTINATION ${UTILS_INCDIR}/io)

#----------------------------------------------------------------
# C) Define the targets for package's library(s)
#----------------------------------------------------------------
TRIBITS_ADD_LIBRARY(
  utils
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
