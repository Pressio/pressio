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
  )

APPEND_SET(HEADERS
  integrators/ode_integrate_n_steps.hpp
  )

APPEND_SET(HEADERS
  step_methods/base/ode_explicit_stepper_base.hpp
  step_methods/base/ode_implicit_stepper_base.hpp
  step_methods/ode_stepper_traits.hpp
)
APPEND_SET(HEADERS
  step_methods/ode_explicit_euler_stepper.hpp
  step_methods/ode_implicit_euler_stepper.hpp
  step_methods/impl/ode_euler_implicit_impl.hpp
)
APPEND_SET(HEADERS
  step_methods/policies/ode_explicit_euler_policy_base.hpp
  step_methods/policies/ode_explicit_euler_standard_policy.hpp
  step_methods/policies/ode_implicit_euler_policy_base.hpp
  step_methods/policies/ode_implicit_euler_standard_policy.hpp
)


# this is temporary until tribits supports header only
APPEND_SET(SOURCES
  src/ode_forward_declarations.cc
  src/ode_void_collector.cc
  integrators/src/ode_integrate_n_steps.cc
  #
  step_methods/src/ode_explicit_stepper_base.cc
  step_methods/src/ode_implicit_stepper_base.cc
  step_methods/src/ode_euler_implicit_impl.cc
  step_methods/src/ode_explicit_euler_stepper.cc
  step_methods/src/ode_implicit_euler_stepper.cc
  step_methods/src/ode_stepper_traits.cc
  #
  step_methods/src/ode_explicit_euler_policy_base.cc
  step_methods/src/ode_explicit_euler_standard_policy.cc
  step_methods/src/ode_implicit_euler_policy_base.cc
  step_methods/src/ode_implicit_euler_standard_policy.cc
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
