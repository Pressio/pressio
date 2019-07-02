
TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES( KOKKOS
  # REQUIRED_HEADERS  SimpleTpl.hpp
  #REQUIRED_LIBS_NAMES libkokkos.a
  #MUST_FIND_ALL_HEADERS
  )

  # we need to strip kokkos from trilinos but for now leave it
  # later on we would like to possibly build wihtout
  # trilinos if needed but kokkos. So we need to have kokkos libs
  # separate from trilinos libraries
  #KOKKOS         "cmake/tpls/"  PT
