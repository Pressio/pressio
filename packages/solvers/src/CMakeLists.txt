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
  #
  solvers_forward_declarations.hpp
  linear/solvers_linear_base.hpp
  linear/solvers_linear_eigen.hpp
  #
)
###########################
## this experimental stuff
###########################
APPEND_SET(HEADERS
  nonlinear/experimental/solvers_nonlinear_newton_raphson.hpp
)


APPEND_SET(SOURCES
  # least_squares/src/algebra_leastsquares.cc
  # src/algebra_newtonraphson.cc
  src/solvers_forward_declarations.cc
  linear/src/solvers_linear_base.cc
  linear/src/solvers_linear_eigen.cc
)
###########################
## this experimental stuff
###########################
APPEND_SET(SOURCES
  nonlinear/experimental/solvers_nonlinear_newton_raphson.cc
)


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
