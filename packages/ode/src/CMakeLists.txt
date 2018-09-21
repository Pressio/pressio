
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------
#To include "ode_config.h"
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ode)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

FILE(GLOB ode_general_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/*.hpp")
FILE(GLOB ode_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/meta/*.hpp")
FILE(GLOB ode_integrators_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/integrators/*.hpp")

#-------------------------
# POLICIES
#------------------------- 
#policies base
FILE(GLOB ode_policies_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/policies/base/*.hpp")
#policies meta
FILE(GLOB ode_policies_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/policies/meta/*.hpp")
#policies standard
FILE(GLOB ode_policies_standard_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/policies/standard/*.hpp")
# custom
FILE(GLOB ode_custom_policies_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/policies/custom/*.hpp")

#-------------------------
# STEPPERS
#-------------------------
FILE(GLOB ode_steppers_common_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
   "${CMAKE_CURRENT_SOURCE_DIR}/ode/steppers/*.hpp")
  
#-------------------------
# EXPLICIT STEPPERS
#-------------------------  
FILE(GLOB ode_explicit_steppers_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/steppers/explicit_steppers/base/*.hpp")
FILE(GLOB ode_explicit_steppers_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/steppers/explicit_steppers/*.hpp")
FILE(GLOB ode_explicit_steppers_impl_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/steppers/explicit_steppers/impl/*.hpp")

#-------------------------
# IMPLICIT STEPPERS
#-------------------------  
FILE(GLOB ode_implicit_steppers_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/steppers/implicit_steppers/base/*.hpp")
FILE(GLOB ode_implicit_steppers_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/steppers/implicit_steppers/*.hpp")
FILE(GLOB ode_implicit_steppers_impl_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/steppers/implicit_steppers/impl/*.hpp")


# SET(HEADERS
#   ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
#   ${ode_meta_HEADERS}
#   ${ode_integrators_HEADERS}
#   #
#   ${ode_policies_base_HEADERS} 
#   ${ode_policies_meta_HEADERS}
#   ${ode_policies_standard_HEADERS}
#   ${ode_custom_policies_HEADERS}
#   #
#   ${ode_steppers_common_HEADERS}
#   ${ode_explicit_steppers_base_HEADERS} 
#   ${ode_explicit_steppers_impl_HEADERS} 
#   ${ode_explicit_steppers_HEADERS} 
#   ${ode_implicit_steppers_base_HEADERS}
#   ${ode_implicit_steppers_HEADERS} 
#   ${ode_implicit_steppers_impl_HEADERS}
# )


#######################
# all sources
#######################

FILE(GLOB ode_general_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/src/*.cc")
FILE(GLOB ode_meta_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/meta/src/*.cc")
FILE(GLOB ode_integrators_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/integrators/src/*.cc")

#-------------------------
# POLICIES
#-------------------------   
#policies  base
FILE(GLOB ode_policies_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/policies/base/src/*.cc")
#policies meta
FILE(GLOB ode_policies_meta_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/policies/meta/src/*.hpp")
#policies standard
FILE(GLOB ode_policies_standard_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/policies/standard/src/*.hpp")
# custom
FILE(GLOB ode_custom_policies_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/policies/custom/src/*.cc")

#-------------------------
# STEPPERS
#-------------------------  
FILE(GLOB ode_steppers_common_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
   "${CMAKE_CURRENT_SOURCE_DIR}/ode/steppers/src/*.cc")
  
#-------------------------
# EXPLICIT STEPPERS
#-------------------------  
FILE(GLOB ode_explicit_steppers_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/steppers/explicit_steppers/base/src/*.cc")
FILE(GLOB ode_explicit_steppers_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/steppers/explicit_steppers/src/*.cc")
FILE(GLOB ode_explicit_steppers_impl_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/steppers/explicit_steppers/impl/src/*.cc")

#-------------------------
# IMPLICIT STEPPERS
#-------------------------  
FILE(GLOB ode_implicit_steppers_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/steppers/implicit_steppers/base/src/*.cc")
FILE(GLOB ode_implicit_steppers_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/steppers/implicit_steppers/src/*.cc")
FILE(GLOB ode_implicit_steppers_impl_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/ode/steppers/implicit_steppers/impl/src/*.cc")


SET(SOURCES
  ${ode_general_SRCS}
  ${ode_meta_SRCS}
  ${ode_integrators_SRCS}
  #
  ${ode_policies_base_SRCS} 
  ${ode_policies_meta_SRCS}
  ${ode_policies_standard_SRCS}
  ${ode_custom_policies_SRCS}
  #
  ${ode_steppers_common_SRCS}
  ${ode_explicit_steppers_base_SRCS} 
  ${ode_explicit_steppers_impl_SRCS} 
  ${ode_explicit_steppers_SRCS} 
  ${ode_implicit_steppers_base_SRCS} 
  ${ode_implicit_steppers_SRCS} 
  ${ode_implicit_steppers_impl_SRCS}
)


SET(ODE_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

#--------------------------------------
SET(LUMPED_HEADERS
  ODE_ALL
  ODE_BASIC
  ODE_INTEGRATORS
  ODE_POLICIES
  ODE_STEPPERS  
)
INSTALL(FILES ${LUMPED_HEADERS} DESTINATION ${ODE_INCDIR})


#--------------------------------------
# install the actual headers
#--------------------------------------
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h DESTINATION ${ODE_INCDIR}/ode)

INSTALL(FILES ${ode_general_HEADERS} DESTINATION ${ODE_INCDIR}/ode)
INSTALL(FILES ${ode_meta_HEADERS} DESTINATION ${ODE_INCDIR}/ode/meta)
INSTALL(FILES ${ode_integrators_HEADERS} DESTINATION ${ODE_INCDIR}/ode/integrators)
#
INSTALL(FILES ${ode_policies_base_HEADERS} DESTINATION ${ODE_INCDIR}/ode/policies/base)
INSTALL(FILES ${ode_policies_meta_HEADERS} DESTINATION ${ODE_INCDIR}/odepolicies/meta)
INSTALL(FILES ${ode_policies_standard_HEADERS} DESTINATION ${ODE_INCDIR}/ode/policies/standard)
INSTALL(FILES ${ode_custom_policies_HEADERS} DESTINATION ${ODE_INCDIR}/ode/policies/custom)
#
INSTALL(FILES ${ode_steppers_common_HEADERS} DESTINATION ${ODE_INCDIR}/ode/steppers)

INSTALL(FILES ${ode_explicit_steppers_base_HEADERS} DESTINATION ${ODE_INCDIR}/ode/steppers/explicit_steppers/base)
INSTALL(FILES ${ode_explicit_steppers_impl_HEADERS} DESTINATION ${ODE_INCDIR}/ode/steppers/explicit_steppers/impl)
INSTALL(FILES ${ode_explicit_steppers_HEADERS} DESTINATION ${ODE_INCDIR}/ode/steppers/explicit_steppers)

INSTALL(FILES ${ode_implicit_steppers_base_HEADERS} DESTINATION ${ODE_INCDIR}/ode/steppers/implicit_steppers/base)
INSTALL(FILES ${ode_implicit_steppers_HEADERS} DESTINATION ${ODE_INCDIR}/ode/steppers/implicit_steppers)
INSTALL(FILES ${ode_implicit_steppers_impl_HEADERS} DESTINATION ${ODE_INCDIR}/ode/steppers/implicit_steppers/impl)


TRIBITS_ADD_LIBRARY(
  ode
  SOURCES ${SOURCES}
  DEPLIBS
  )


# #----------------------------------------------------------------
# # C) Define the targets for package's library(s)
# #----------------------------------------------------------------
# TRIBITS_ADD_LIBRARY(
#   ode
#   HEADERS ${HEADERS}
#   SOURCES ${SOURCES}
# )

# # D) Define an executable that gets installed
# TRIBITS_ADD_EXECUTABLE(commomBase-helloworld  NOEXEPREFIX  NOEXESUFFIX
#   SOURCES  SimpleCxx_HelloWorld_Prog.cpp
#   INSTALLABLE
#   )

