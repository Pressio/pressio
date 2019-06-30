
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}) # To include "ode_config.h"

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------

#######################
# HEADERS
#######################

FILE(GLOB ode_general_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
FILE(GLOB ode_integrators_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/integrators/*.hpp")

FILE(GLOB ode_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/meta/*.hpp")

# EXPLICIT
FILE(GLOB ode_exp_policies_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/explicit/policies/*.hpp")
FILE(GLOB ode_exp_steppers_impl_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
   "${CMAKE_CURRENT_SOURCE_DIR}/explicit/steppers/impl/*.hpp")
FILE(GLOB ode_exp_steppers_common_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
   "${CMAKE_CURRENT_SOURCE_DIR}/explicit/steppers/*.hpp")

# IMPLICIT
FILE(GLOB ode_imp_comm_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/implicit/*.hpp")
FILE(GLOB ode_imp_policies_base_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/implicit/policies/base/*.hpp")
FILE(GLOB ode_imp_policies_meta_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/implicit/policies/meta/*.hpp")
FILE(GLOB ode_imp_policies_std_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/implicit/policies/standard/*.hpp")
FILE(GLOB ode_imp_steppers_common_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
   "${CMAKE_CURRENT_SOURCE_DIR}/implicit/steppers/*.hpp")

# all sources
FILE(GLOB ode_general_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.cc")
SET(SOURCES ${ode_general_SRCS})

#--------------------------------------
SET(LUMPED_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ODE_ALL
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ODE_BASIC
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ODE_INTEGRATORS
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ODE_EXPLICIT
  ${CMAKE_CURRENT_SOURCE_DIR}/../../ODE_IMPLICIT
)
INSTALL(FILES ${LUMPED_HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

#--------------------------------------
# install the actual headers
#--------------------------------------
SET(ODE_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/ode/src)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h DESTINATION ${ODE_INCDIR})

INSTALL(FILES ${ode_general_HEADERS} DESTINATION ${ODE_INCDIR})
INSTALL(FILES ${ode_integrators_HEADERS} DESTINATION ${ODE_INCDIR}/integrators)
INSTALL(FILES ${ode_meta_HEADERS} DESTINATION ${ODE_INCDIR}/meta)

INSTALL(FILES ${ode_exp_policies_HEADERS} DESTINATION ${ODE_INCDIR}/explicit/policies)
INSTALL(FILES ${ode_exp_steppers_impl_HEADERS} DESTINATION ${ODE_INCDIR}/explicit/steppers/impl)
INSTALL(FILES ${ode_exp_steppers_common_HEADERS} DESTINATION ${ODE_INCDIR}/explicit/steppers)

INSTALL(FILES ${ode_imp_comm_HEADERS} DESTINATION ${ODE_INCDIR}/implicit)
INSTALL(FILES ${ode_imp_policies_base_HEADERS} DESTINATION ${ODE_INCDIR}/implicit/policies/base)
INSTALL(FILES ${ode_imp_policies_meta_HEADERS} DESTINATION ${ODE_INCDIR}/implicit/policies/meta)
INSTALL(FILES ${ode_imp_policies_std_HEADERS} DESTINATION ${ODE_INCDIR}/implicit/policies/standard)
INSTALL(FILES ${ode_imp_steppers_common_HEADERS} DESTINATION ${ODE_INCDIR}/implicit/steppers)

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
