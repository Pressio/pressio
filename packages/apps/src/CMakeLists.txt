
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}) #To include "name_config.h"

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------

FILE(GLOB apps_gen_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/*.hpp")

FILE(GLOB apps_burgers1d_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/burgers1d/*.hpp")

FILE(GLOB apps_steadyLinearAdvDiff1d_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/steady_linear_adv_diff1d/*.hpp")

FILE(GLOB apps_steadyLinearAdvDiff2d_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/steady_linear_adv_diff2d/*.hpp")

FILE(GLOB apps_unsteadyLinearAdvDiff1d_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/unsteady_linear_adv_diff1d/*.hpp")

FILE(GLOB apps_unsteadyNonLinearAdvDiffReaction2d_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/unsteady_nonlinear_adv_diff_reaction2d/*.hpp")

FILE(GLOB apps_unsteadyNonLinearAdvDiffReactionFlame2d_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/unsteady_nonlinear_adv_diff_reaction_flame2d/*.hpp")

# sources
FILE(GLOB apps_burgers1d_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/burgers1d/*.cc")

FILE(GLOB apps_steadyLinearAdvDiff1d_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/steady_linear_adv_diff1d/*.cc")

FILE(GLOB apps_steadyLinearAdvDiff2d_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/steady_linear_adv_diff2d/*.cc")

FILE(GLOB apps_unsteadyLinearAdvDiff1d_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/unsteady_linear_adv_diff1d/*cc")

FILE(GLOB apps_unsteadyNonLinearAdvDiffReaction2d_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/unsteady_nonlinear_adv_diff_reaction2d/*.cc")

FILE(GLOB apps_unsteadyNonLinearAdvDiffReactionFlame2d_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
  "${CMAKE_CURRENT_SOURCE_DIR}/unsteady_nonlinear_adv_diff_reaction_flame2d/*.cc")

SET(SOURCES
  ${apps_burgers1d_SRCS}
  ${apps_steadyLinearAdvDiff1d_SRCS}
  ${apps_steadyLinearAdvDiff2d_SRCS}
  ${apps_unsteadyLinearAdvDiff1d_SRCS}
  ${apps_unsteadyNonLinearAdvDiffReaction2d_SRCS}
  ${apps_unsteadyNonLinearAdvDiffReactionFlame2d_SRCS}
)

#--------------------------------------
SET(LUMPED_HEADERS
  ${CMAKE_CURRENT_SOURCE_DIR}/../../APPS_BURGERS1D
  ${CMAKE_CURRENT_SOURCE_DIR}/../../APPS_STEADYLINADVDIFF1D
  ${CMAKE_CURRENT_SOURCE_DIR}/../../APPS_STEADYLINADVDIFF2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../../APPS_UNSTEADYLINADVDIFF1D
  ${CMAKE_CURRENT_SOURCE_DIR}/../../APPS_UNSTEADYNONLINADVDIFFREACTION2D
  ${CMAKE_CURRENT_SOURCE_DIR}/../../APPS_UNSTEADYNONLINADVDIFFREACTIONFLAME2D
)
INSTALL(FILES ${LUMPED_HEADERS} DESTINATION ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

#--------------------------------------
# install the actual headers
#--------------------------------------
SET(APPS_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR}/apps/src)

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h DESTINATION ${APPS_INCDIR})
INSTALL(FILES ${apps_gen_HEADERS} DESTINATION ${APPS_INCDIR})
INSTALL(FILES ${apps_burgers1d_HEADERS} DESTINATION ${APPS_INCDIR}/burgers1d)
INSTALL(FILES ${apps_steadyLinearAdvDiff1d_HEADERS} DESTINATION ${APPS_INCDIR}/steady_linear_adv_diff1d)
INSTALL(FILES ${apps_steadyLinearAdvDiff2d_HEADERS} DESTINATION ${APPS_INCDIR}/steady_linear_adv_diff2d)
INSTALL(FILES ${apps_unsteadyLinearAdvDiff1d_HEADERS} DESTINATION ${APPS_INCDIR}/unsteady_linear_adv_diff1d)
INSTALL(FILES ${apps_unsteadyNonLinearAdvDiffReaction2d_HEADERS} DESTINATION ${APPS_INCDIR}/unsteady_nonlinear_adv_diff_reaction2d)
INSTALL(FILES ${apps_unsteadyNonLinearAdvDiffReactionFlame2d_HEADERS} DESTINATION ${APPS_INCDIR}/unsteady_nonlinear_adv_diff_reactionFlame2d)

#----------------------------------------------------------------
# C) Define the targets for package's library(s)
#----------------------------------------------------------------
TRIBITS_ADD_LIBRARY(
  apps
  #   HEADERS ${HEADERS}
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
