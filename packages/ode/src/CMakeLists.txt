#
# A) Package-specific configuration options
#
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

# B) Define the header and source files (and include directories)
SET(HEADERS "")
SET(SOURCES "")

#To include "ode_config.h" (because it lives in the binary folder)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
# # To include the configured header above
# INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

#######################
# include all headers
#######################

FILE(GLOB ode_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/meta/*.hpp")

# integrators
FILE(GLOB ode_integrators_impl_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/integrators/impl/*.hpp")
FILE(GLOB ode_integrators_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/integrators/*.hpp")

# policies common
FILE(GLOB ode_policies_common_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
     "${CMAKE_CURRENT_SOURCE_DIR}/policies/common/*.hpp")

#policies base
FILE(GLOB ode_policies_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/policies/base/*.hpp")

#-------------------------
# EXPLICIT POLICIES
#-------------------------  
# euler
FILE(GLOB ode_explicit_policies_euler_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/policies/explicit_policies/euler/*.hpp")
# rk
FILE(GLOB ode_explicit_policies_rk_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/policies/explicit_policies/runge_kutta/*.hpp")
#-------------------------
# EXPLICIT STEPPERS
#-------------------------  
FILE(GLOB ode_explicit_steppers_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/explicit_steppers/base/*.hpp")
FILE(GLOB ode_explicit_steppers_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/explicit_steppers/*.hpp")
FILE(GLOB ode_explicit_steppers_impl_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/explicit_steppers/impl/*.hpp")

#-------------------------
# IMPLICIT POLICIES
#-------------------------  
# euler
FILE(GLOB ode_implicit_policies_euler_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/policies/implicit_policies/euler/*.hpp")
#-------------------------
# IMPLICIT STEPPERS
#-------------------------  
FILE(GLOB ode_implicit_steppers_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/implicit_steppers/base/*.hpp")
FILE(GLOB ode_implicit_steppers_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/implicit_steppers/*.hpp")
FILE(GLOB ode_implicit_steppers_impl_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/implicit_steppers/impl/*.hpp")

#-------------------------
# CUSTOM POLICIES
#------------------------- 
FILE(GLOB ode_custom_policies_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/policies/custom/*.hpp")
  

 
#######################
#######################
# include all sources
#######################
#######################
FILE(GLOB ode_meta_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/meta/src/*.cc")

# integrators
FILE(GLOB ode_integrators_impl_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/integrators/impl/src/*.cc")
FILE(GLOB ode_integrators_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/integrators/src/*.cc")

# policies common
FILE(GLOB ode_policies_common_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
     "${CMAKE_CURRENT_SOURCE_DIR}/policies/common/src/*.cc")
#policies  base
FILE(GLOB ode_policies_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/policies/base/src/*.cc")

#-------------------------
# EXPLICIT POLICIES
#-------------------------  
#policies euler
FILE(GLOB ode_explicit_policies_euler_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/policies/explicit_policies/euler/src/*.cc")
#policies runge_kutta
FILE(GLOB ode_explicit_policies_rk_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/policies/explicit_policies/runge_kutta/src/*.cc")

#-------------------------
# EXPLICIT STEPPERS
#-------------------------  
FILE(GLOB ode_explicit_steppers_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/explicit_steppers/base/src/*.cc")
FILE(GLOB ode_explicit_steppers_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/explicit_steppers/src/*.cc")
FILE(GLOB ode_explicit_steppers_impl_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/explicit_steppers/impl/src/*.cc")

#-------------------------
# IMPLICIT POLICIES
#-------------------------  
#policies euler
FILE(GLOB ode_implicit_policies_euler_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/policies/implicit_policies/euler/src/*.cc")

#-------------------------
# IMPLICIT STEPPERS
#-------------------------  
FILE(GLOB ode_implicit_steppers_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/implicit_steppers/base/src/*.cc")
FILE(GLOB ode_implicit_steppers_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/implicit_steppers/src/*.cc")
FILE(GLOB ode_implicit_steppers_impl_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/implicit_steppers/impl/src/*.cc")

#-------------------------
# CUSTOM POLICIES
#------------------------- 
FILE(GLOB ode_custom_policies_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/policies/custom/src/*.cc")
  

  

APPEND_SET(HEADERS
  ${ode_meta_HEADERS}
  ${ode_integrators_impl_HEADERS}
  ${ode_integrators_HEADERS}
  #
  ${ode_policies_common_HEADERS}
  ${ode_policies_base_HEADERS} 
  ${ode_explicit_policies_euler_HEADERS}
  ${ode_explicit_policies_rk_HEADERS}
  ${ode_implicit_policies_euler_HEADERS}
  ${ode_custom_policies_HEADERS}
  #
  ${ode_explicit_steppers_base_HEADERS} 
  ${ode_explicit_steppers_impl_HEADERS} 
  ${ode_explicit_steppers_HEADERS} 
  ${ode_implicit_steppers_base_HEADERS}
  ${ode_implicit_steppers_HEADERS} 
  ${ode_implicit_steppers_impl_HEADERS}
  ${ode_custom_policies_SRCS}
)

APPEND_SET(SOURCES
  ${ode_meta_SRCS}
  ${ode_integrators_impl_SRCS}
  ${ode_integrators_SRCS}
  ${ode_policies_common_SRCS}
  ${ode_policies_base_SRCS}
  ${ode_explicit_policies_euler_SRCS}
  ${ode_explicit_policies_rk_SRCS}
  ${ode_implicit_policies_euler_SRCS}
  #
  ${ode_explicit_steppers_base_SRCS} 
  ${ode_explicit_steppers_impl_SRCS} 
  ${ode_explicit_steppers_SRCS} 
  ${ode_implicit_steppers_base_SRCS} 
  ${ode_implicit_steppers_SRCS} 
  ${ode_implicit_steppers_impl_SRCS}
)

# C) Define the targets for package's library(s)
TRIBITS_ADD_LIBRARY(
  ode
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )
# # this seems to work for header only library
# TRIBITS_INSTALL_HEADERS(
#   HEADERS ${HEADERS}
# )

# # D) Define an executable that gets installed
# TRIBITS_ADD_EXECUTABLE(commomBase-helloworld  NOEXEPREFIX  NOEXESUFFIX
#   SOURCES  SimpleCxx_HelloWorld_Prog.cpp
#   INSTALLABLE
#   )














# #-------------------------
# # IMPLICIT POLICIES
# #-------------------------
# #policies implicit base
# FILE(GLOB ode_implicit_policies_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/policies/implicit_policies/base/*.hpp")
# #policies euler
# FILE(GLOB ode_implicit_policies_euler_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/policies/implicit_policies/euler/*.hpp")
# #policies bdf2
# FILE(GLOB ode_implicit_policies_bdf2_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/policies/implicit_policies/bdf2/*.hpp")
# #policies adams-moulton1
# FILE(GLOB ode_implicit_policies_adams_moulton1_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/policies/implicit_policies/adams_moulton/*.hpp")

# #-------------------------
# # IMPLICIT STEPPERS
# #-------------------------  
# FILE(GLOB ode_implicit_steppers_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/steppers/implicit_steppers/base/*.hpp")
# FILE(GLOB ode_implicit_steppers_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/steppers/implicit_steppers/*.hpp")
# FILE(GLOB ode_implicit_steppers_impl_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/steppers/implicit_steppers/impl/*.hpp")

  


  # #-------------------------
# # IMPLICIT STEPPERS
# #-------------------------  
# FILE(GLOB ode_implicit_steppers_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/steppers/implicit_steppers/base/src/*.cc")
# FILE(GLOB ode_implicit_steppers_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/steppers/implicit_steppers/src/*.cc")
# FILE(GLOB ode_implicit_steppers_impl_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/steppers/implicit_steppers/impl/src/*.cc")



# #-------------------------
# # IMPLICIT POLICIES
# #-------------------------
# #policies implicit base
# FILE(GLOB ode_implicit_policies_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/policies/implicit_policies/base/src/*.cc")
# #policies euler
# FILE(GLOB ode_implicit_policies_euler_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/policies/implicit_policies/euler/src/*.cc")
# #policies bdf2
# FILE(GLOB ode_implicit_policies_bdf2_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/policies/implicit_policies/bdf2/src/*.cc")
# #policies adams-moulton1
# FILE(GLOB ode_implicit_policies_adams_moulton1_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#    "${CMAKE_CURRENT_SOURCE_DIR}/policies/implicit_policies/adams_moulton/src/*.cc")

