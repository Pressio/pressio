
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------
# To include solvers_config.h
# INCLUDE_DIRECTORIES()
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/solvers)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

# ###############
# # HEADERS
# ###############
# FILE(GLOB solvers_general_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/solvers/*.hpp")
# FILE(GLOB solvers_exp_HEADERS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
#     "${CMAKE_CURRENT_SOURCE_DIR}/solvers/experimental/*.hpp")

SET(GEN_HEADERS
solvers/solvers_ConfigDefs.hpp
solvers/solvers_forward_declarations.hpp
)

SET(EXP_HEADERS
solvers/experimental/solvers_nonlinear_factory.hpp
solvers/experimental/solvers_l2_vector_norm.hpp
solvers/experimental/solvers_linear_base.hpp
solvers/experimental/solvers_linear_direct.hpp
solvers/experimental/solvers_linear_factory.hpp
solvers/experimental/solvers_linear_iterative_policies_trilinos.hpp
solvers/experimental/solvers_linear_iterative.hpp
solvers/experimental/solvers_linear_traits.hpp
solvers/experimental/solvers_meta_static_checks.hpp
solvers/experimental/solvers_nonlinear_base.hpp
solvers/experimental/solvers_nonlinear_iterative.hpp
solvers/experimental/solvers_nonlinear_traits.hpp
solvers/experimental/solvers_norms_fwd.hpp
solvers/experimental/solvers_policy_linear_dense_eigen.hpp
solvers/experimental/solvers_policy_linear_iterative_eigen.hpp
solvers/experimental/solvers_policy_linear_iterative_trilinos.hpp
solvers/experimental/solvers_policy_nonlinear_iterative.hpp
solvers/experimental/system_traits.hpp
)


#################
#### SOURCES
#################
FILE(GLOB solvers_general_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/solvers/src/*.cc")
FILE(GLOB solvers_exp_SRCS RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_SOURCE_DIR}/solvers/experimental/src/*.cc")

SET(SOURCES 
  ${solvers_general_SRCS}
  ${solvers_exp_SRCS}
)


SET(SOLVERS_INCDIR ${CMAKE_INSTALL_PREFIX}/${${PROJECT_NAME}_INSTALL_INCLUDE_DIR})

#--------------------------------------
SET(LUMPED_HEADERS
  SOLVERS_ALL
  SOLVERS_BASIC
  SOLVERS_LINEAR
  SOLVERS_NONLINEAR  
)
INSTALL(FILES ${LUMPED_HEADERS} DESTINATION ${SOLVERS_INCDIR})

#--------------------------------------
# install the actual core headers
#--------------------------------------
INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h DESTINATION ${SOLVERS_INCDIR}/solvers)

INSTALL(FILES ${GEN_HEADERS} DESTINATION ${SOLVERS_INCDIR}/solvers)
INSTALL(FILES ${EXP_HEADERS} DESTINATION ${SOLVERS_INCDIR}/solvers/experimental)


#----------------------------------------------------------------
# C) Define the targets for package's library(s)
#----------------------------------------------------------------
TRIBITS_ADD_LIBRARY(
  solvers
  SOURCES ${SOURCES}
  DEPLIBS
  )
# TRIBITS_ADD_LIBRARY(
#   solvers
#   HEADERS ${HEADERS}
#   SOURCES ${SOURCES}
#   )

# #
# # D) Define an executable that gets installed
# #
# TRIBITS_ADD_EXECUTABLE(commomBase-helloworld  NOEXEPREFIX  NOEXESUFFIX
#   SOURCES  SimpleCxx_HelloWorld_Prog.cpp
#   INSTALLABLE
#   )
