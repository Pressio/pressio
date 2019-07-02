
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}) #To include "mpl_config.h"

# Set the base directory in order to find pressio_version.h
ASSERT_DEFINED(PROJECT_SOURCE_DIR)
INCLUDE_DIRECTORIES(${PROJECT_NAME}_BINARY_DIR})


#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------
SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  ${CMAKE_CURRENT_SOURCE_DIR}/optimization_ConfigDefs.hpp
  ${CMAKE_CURRENT_SOURCE_DIR}/optimization_version.hpp
)

SET(SOURCES
  ${CMAKE_CURRENT_SOURCE_DIR}/optimization_version.cc
)

#----------------------------------------------------------------
# C) Define the targets for package's library(s)
#----------------------------------------------------------------
TRIBITS_ADD_LIBRARY(
  optimization
  # HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

# #----------------------------------------------------------------
# # D) Define an executable that gets installed
# #----------------------------------------------------------------
# TRIBITS_ADD_EXECUTABLE(commomBase-helloworld  NOEXEPREFIX  NOEXESUFFIX
#   SOURCES  SimpleCxx_HelloWorld_Prog.cpp
#   INSTALLABLE
#   )
