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
  ode_ConfigDefs.hpp
  ode_forward_declarations.hpp
  ode_void_collector.hpp
  ode_stepper_traits.hpp
  ode_explicit_stepper_base.hpp
  ode_euler_stepper.hpp
  ode_rk4_stepper.hpp
  ode_integrate_n_steps.hpp
)

APPEND_SET(SOURCES
  src/ode_forward_declarations.cc
  src/ode_void_collector.cc
  src/ode_stepper_traits.cc  
  src/ode_explicit_stepper_base.cc
  src/ode_euler_stepper.cc
  src/ode_rk4_stepper.cc
  src/ode_integrate_n_steps.cc
  )

#
# C) Define the targets for package's library(s)
#

# # this seems to work for header only library
# TRIBITS_INSTALL_HEADERS(
#   HEADERS ${HEADERS}
# )
# this is the original in case we have sources too for the library
TRIBITS_ADD_LIBRARY(
  ode
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
