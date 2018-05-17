TRIBITS_REPOSITORY_DEFINE_TPLS(
  MPI             "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/FindTPLMPI.cmake"   PT
  TRILINOS        "cmake/tpls/"  PT
  EIGEN           "cmake/tpls/"  PT
  )

# PT means required 
# ST optional 
