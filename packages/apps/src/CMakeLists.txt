
# A) Package-specific configuration options
TRIBITS_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#----------------------------------------------------------------
# B) Define the header and source files (and include directories)
#----------------------------------------------------------------
# To include "apps_config.h"
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

SET(HEADERS
  apps_ConfigDefs.hpp
  apps_burgers1d_eigen.hpp
  apps_burgers1d_epetra.hpp
  #
  apps_constants.hpp
  apps_mpi_square_grid.hpp
  apps_mpi_square_grid.hpp
  apps_square_periodic_cc_mesh.hpp
  apps_taylor_green.hpp
  apps_shear_layer.hpp
)

SET(SOURCES
  src/apps_burgers1d_eigen.cc
  src/apps_burgers1d_epetra.cc
  src/apps_constants.cc
  src/apps_mpi_square_grid.cc
  src/apps_square_periodic_cc_mesh.cc
  src/apps_taylor_green.cc
  src/apps_shear_layer.cc
)

#----------------------------------------------------------------
# C) Define the targets for package's library(s)
#----------------------------------------------------------------
TRIBITS_ADD_LIBRARY(
  apps
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
