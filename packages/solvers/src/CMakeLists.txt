
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------
# To include solvers_config.h
# INCLUDE_DIRECTORIES()
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  solvers_ConfigDefs.hpp
  solvers_forward_declarations.hpp
)
SET(SOURCES
  src/solvers_forward_declarations.cc
)

# experimental/hacked stuff 
APPEND_SET(HEADERS
  experimental/solvers_linear_iterative.hpp
  experimental/solvers_linear_factory.hpp
  )

APPEND_SET(SOURCES
  experimental/src/solvers_linear_iterative.cc
  experimental/src/solvers_linear_factory.cc
)


#----------------------------------------------------------------
# C) Define the targets for package's library(s)
#----------------------------------------------------------------
TRIBITS_ADD_LIBRARY(
  solvers
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
