
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}) #To include "core_config.h"

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------

# ###############
# # HEADERS
# ###############
# FILE(GLOB solvers_general_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/solvers/*.hpp")
# FILE(GLOB solvers_exp_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/solvers/experimental/*.hpp")

SET(GEN_HEADERS
solvers_ConfigDefs.hpp
solvers_forward_declarations.hpp
)

SET(EXP_HEADERS
experimental/solvers_l2_vector_norm.hpp
experimental/solvers_linear_base.hpp
experimental/solvers_linear_direct.hpp
experimental/solvers_linear_factory.hpp
experimental/solvers_linear_iterative_policies_trilinos.hpp
experimental/solvers_linear_iterative.hpp
experimental/solvers_linear_traits.hpp
experimental/solvers_meta_static_checks.hpp
experimental/solvers_nonlinear_base.hpp
experimental/solvers_nonlinear_factory.hpp
experimental/solvers_nonlinear_iterative.hpp
experimental/solvers_nonlinear_traits.hpp
experimental/solvers_norms_fwd.hpp
experimental/solvers_policy_linear_dense_eigen.hpp
experimental/solvers_policy_linear_iterative_eigen.hpp
experimental/solvers_policy_linear_iterative_trilinos.hpp
experimental/solvers_policy_nonlinear_iterative.hpp
experimental/solvers_policy_nonlinear_leastsquare_iterative.hpp
experimental/solvers_system_traits.hpp
)


#################
#### SOURCES
#################
FILE(GLOB solvers_general_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc")
FILE(GLOB solvers_exp_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/experimental/src/*.cc")

SET(SOURCES 
  ${solvers_general_SRCS}
  ${solvers_exp_SRCS}
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
INSTALL(FILES ${EXP_HEADERS} DESTINATION ${SOLVERS_INCDIR}/experimental)


#----------------------------------------------------------------
# C) Define the targets for package's library(s)
#----------------------------------------------------------------
TRIBITS_ADD_LIBRARY(
  solvers
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
