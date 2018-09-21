
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------
# To include optimization_config.h
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  optimization_ConfigDefs.hpp
  dummy.hpp
  #
  #least_squares/algebra_leastsquares.hpp
  #algebra_newtonraphson.hpp
)

SET(SOURCES
  dummy.cc
  #least_squares/src/algebra_leastsquares.cc
  #src/algebra_newtonraphson.cc
)

#----------------------------------------------------------------
# C) Define the targets for package's library(s)
#----------------------------------------------------------------
TRIBITS_ADD_LIBRARY(
  optimization
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

# #----------------------------------------------------------------
# # D) Define an executable that gets installed
# #----------------------------------------------------------------
# TRIBITS_ADD_EXECUTABLE(commomBase-helloworld  NOEXEPREFIX  NOEXESUFFIX
#   SOURCES  SimpleCxx_HelloWorld_Prog.cpp
#   INSTALLABLE
#   )
