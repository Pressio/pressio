TRIBITS_REPOSITORY_DEFINE_TPLS(
  MPI             "cmake/tpls/"  PT
  GTEST           "cmake/tpls/"  PT
  TRILINOS        "cmake/tpls/"  PT  
  EIGEN           "cmake/tpls/"  PT
  BLAZE           "cmake/tpls/"  PT

  # we need to strip kokkos from trilinos but for now leave it
  # later on we would like to possibly build rompp wihtout
  # trilinos if needed but kokkos. So we need to have kokkos libs
  # separate from trilinos libraries
  #KOKKOS         "cmake/tpls/"  PT
  )

# PT means required 
# ST optional 
#MPI            "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/FindTPLMPI.cmake"   PT
