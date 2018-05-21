#
# A) Package-specific configuration options
#

TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

# To include the configured header above
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  svd_ConfigDefs.hpp
  svd_forward_declarations.hpp
  #
  svd_solver_eigen.hpp
  svd_solver_generic_base.hpp
  svd_solver_traits.hpp  
)

APPEND_SET(SOURCES
  src/svd_forward_declarations.cc
  src/svd_solver_eigen.cc
  src/svd_solver_generic_base.cc
  src/svd_solver_traits.cc
)

#
# C) Define the targets for package's library(s)
#

# # # this seems to work for header only library
# TRIBITS_INSTALL_HEADERS(
#   HEADERS ${HEADERS}
# )
# this is the original in case we have sources too for the library
TRIBITS_ADD_LIBRARY(
  svd
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
