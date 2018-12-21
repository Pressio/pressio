
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}) #To include "core_config.h"

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------

# ###############
# # HEADERS
# ###############
FILE(GLOB GEN_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")
FILE(GLOB LINEAR_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/linear/*.hpp")

#################
#### SOURCES
#################
FILE(GLOB solvers_general_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc")
FILE(GLOB solvers_linear_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/linear/src/*.cc")

SET(SOURCES
  ${solvers_general_SRCS}
  # ${solvers_exp_SRCS}
)

#--------------------------------------
SET(LUMPED_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/../../SOLVERS_ALL
  ${CMAKE_CURRENT_SOURCE_DIR}/../../SOLVERS_BASIC
  ${CMAKE_CURRENT_SOURCE_DIR}/../../SOLVERS_LINEAR
  ${CMAKE_CURRENT_SOURCE_DIR}/../../SOLVERS_NONLINEAR
)
INSTALL(FILES ${LUMPED_HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})


#--------------------------------------
# install the actual core headers
#--------------------------------------
SET(SOLVERS_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/solvers/src)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h DESTINATION ${SOLVERS_INCDIR})
INSTALL(FILES ${GEN_HEADERS} DESTINATION ${SOLVERS_INCDIR})
INSTALL(FILES ${EXP_HEADERS} DESTINATION ${SOLVERS_INCDIR}/linear)

#----------------------------------------------------------------
# C) Define the targets for package's library(s)
#----------------------------------------------------------------
TRIBITS_ADD_LIBRARY(
  solvers
  SOURCES ${SOURCES}
  DEPLIBS
  )

# #
# # D) Define an executable that gets installed
# #
# TRIBITS_ADD_EXECUTABLE(commomBase-helloworld  NOEXEPREFIX  NOEXESUFFIX
#   SOURCES  SimpleCxx_HelloWorld_Prog.cpp
#   INSTALLABLE
#   )
