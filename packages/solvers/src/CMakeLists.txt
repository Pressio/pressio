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
  solvers_ConfigDefs.hpp
  solvers_forward_declarations.hpp
)
APPEND_SET(SOURCES
  src/solvers_forward_declarations.cc
)


########################################################
##           experimental/hacked stuff 
########################################################
APPEND_SET(HEADERS
  experimental/solvers_linear_base.hpp
  experimental/solvers_linear_eigen.hpp
  experimental/solvers_nonlinear_newton_raphson.hpp
)
APPEND_SET(SOURCES
  experimental/src/solvers_linear_base.cc
  experimental/src/solvers_linear_eigen.cc
  experimental/src/solvers_nonlinear_newton_raphson.cc
)
########################################################

# # # this seems to work for header only library
# TRIBITS_INSTALL_HEADERS(
#   HEADERS ${HEADERS}
# )
# this is the original in case we have sources too for the library
TRIBITS_ADD_LIBRARY(
  solvers
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
