
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}) # To include "ode_config.h"

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------

FILE(GLOB ode_general_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
FILE(GLOB ode_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/meta/*.hpp")
FILE(GLOB ode_integrators_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/integrators/*.hpp")

#-------------------------
# POLICIES
#policies base
FILE(GLOB ode_policies_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/policies/base/*.hpp")
#policies meta
FILE(GLOB ode_policies_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/policies/meta/*.hpp")
#policies standard
FILE(GLOB ode_policies_standard_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/policies/standard/*.hpp")
# custom
FILE(GLOB ode_custom_policies_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/policies/custom/*.hpp")

#-------------------------
# STEPPERS
FILE(GLOB ode_steppers_common_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
   "${CMAKE_CURRENT_SOURCE_DIR}/steppers/*.hpp")
  
#-------------------------
# EXPLICIT STEPPERS
FILE(GLOB ode_explicit_steppers_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/explicit_steppers/base/*.hpp")
FILE(GLOB ode_explicit_steppers_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/explicit_steppers/*.hpp")
FILE(GLOB ode_explicit_steppers_impl_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/explicit_steppers/impl/*.hpp")

#-------------------------
# IMPLICIT STEPPERS
FILE(GLOB ode_implicit_steppers_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/implicit_steppers/base/*.hpp")
FILE(GLOB ode_implicit_steppers_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/implicit_steppers/*.hpp")
FILE(GLOB ode_implicit_steppers_impl_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/implicit_steppers/impl/*.hpp")


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
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc")
FILE(GLOB ode_meta_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/meta/src/*.cc")
FILE(GLOB ode_integrators_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/integrators/src/*.cc")

#-------------------------
# POLICIES
#policies  base
FILE(GLOB ode_policies_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/policies/base/src/*.cc")
#policies meta
FILE(GLOB ode_policies_meta_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/policies/meta/src/*.hpp")
#policies standard
FILE(GLOB ode_policies_standard_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/policies/standard/src/*.hpp")
# custom
FILE(GLOB ode_custom_policies_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/policies/custom/src/*.cc")

#-------------------------
# STEPPERS
FILE(GLOB ode_steppers_common_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
   "${CMAKE_CURRENT_SOURCE_DIR}/steppers/src/*.cc")
  
#-------------------------
# EXPLICIT STEPPERS
FILE(GLOB ode_explicit_steppers_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/explicit_steppers/base/src/*.cc")
FILE(GLOB ode_explicit_steppers_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/explicit_steppers/src/*.cc")
FILE(GLOB ode_explicit_steppers_impl_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/explicit_steppers/impl/src/*.cc")

#-------------------------
# IMPLICIT STEPPERS
FILE(GLOB ode_implicit_steppers_base_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/implicit_steppers/base/src/*.cc")
FILE(GLOB ode_implicit_steppers_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/implicit_steppers/src/*.cc")
FILE(GLOB ode_implicit_steppers_impl_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/steppers/implicit_steppers/impl/src/*.cc")


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


#--------------------------------------
SET(LUMPED_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ODE_ALL
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ODE_BASIC
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ODE_INTEGRATORS
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ODE_POLICIES
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ODE_STEPPERS  
)
INSTALL(FILES ${LUMPED_HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})


#--------------------------------------
# install the actual headers
#--------------------------------------
SET(ODE_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/ode/src)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h DESTINATION ${ODE_INCDIR})

INSTALL(FILES ${ode_general_HEADERS} DESTINATION ${ODE_INCDIR})
INSTALL(FILES ${ode_meta_HEADERS} DESTINATION ${ODE_INCDIR}/meta)
INSTALL(FILES ${ode_integrators_HEADERS} DESTINATION ${ODE_INCDIR}/integrators)
#
INSTALL(FILES ${ode_policies_base_HEADERS} DESTINATION ${ODE_INCDIR}/policies/base)
INSTALL(FILES ${ode_policies_meta_HEADERS} DESTINATION ${ODE_INCDIR}/policies/meta)
INSTALL(FILES ${ode_policies_standard_HEADERS} DESTINATION ${ODE_INCDIR}/policies/standard)
INSTALL(FILES ${ode_custom_policies_HEADERS} DESTINATION ${ODE_INCDIR}/policies/custom)
#
INSTALL(FILES ${ode_steppers_common_HEADERS} DESTINATION ${ODE_INCDIR}/steppers)

INSTALL(FILES ${ode_explicit_steppers_base_HEADERS} DESTINATION ${ODE_INCDIR}/steppers/explicit_steppers/base)
INSTALL(FILES ${ode_explicit_steppers_impl_HEADERS} DESTINATION ${ODE_INCDIR}/steppers/explicit_steppers/impl)
INSTALL(FILES ${ode_explicit_steppers_HEADERS} DESTINATION ${ODE_INCDIR}/steppers/explicit_steppers)

INSTALL(FILES ${ode_implicit_steppers_base_HEADERS} DESTINATION ${ODE_INCDIR}/steppers/implicit_steppers/base)
INSTALL(FILES ${ode_implicit_steppers_HEADERS} DESTINATION ${ODE_INCDIR}/steppers/implicit_steppers)
INSTALL(FILES ${ode_implicit_steppers_impl_HEADERS} DESTINATION ${ODE_INCDIR}/steppers/implicit_steppers/impl)


#----------------------------------------------------------------
# C) Define the targets for package's library(s)
#----------------------------------------------------------------
TRIBITS_ADD_LIBRARY(
  ode
  #   HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  DEPLIBS
  )

# # D) Define an executable that gets installed
# TRIBITS_ADD_EXECUTABLE(commomBase-helloworld  NOEXEPREFIX  NOEXESUFFIX
#   SOURCES  SimpleCxx_HelloWorld_Prog.cpp
#   INSTALLABLE
#   )

