
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------
#To include "svd_config.h"
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/svd)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

###############
# HEADERS
###############
FILE(GLOB svd_general_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/svd/*.hpp")

#################
#### SOURCES
#################
FILE(GLOB svd_general_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/svd/src/*.cc")


SET(SOURCES 
  ${svd_general_SRCS}
)


SET(SVD_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

#--------------------------------------
SET(LUMPED_HEADERS
  SVD_BASIC
)
INSTALL(FILES ${LUMPED_HEADERS} DESTINATION ${SVD_INCDIR})

#--------------------------------------
# install the actual headers
#--------------------------------------
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h DESTINATION ${SVD_INCDIR}/svd)

INSTALL(FILES ${svd_general_HEADERS} DESTINATION ${SVD_INCDIR}/svd)


#----------------------------------------------------------------
# C) Define the targets for package's library(s)
#----------------------------------------------------------------
TRIBITS_ADD_LIBRARY(
  svd
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
