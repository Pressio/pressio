TRIBITS_REPOSITORY_DEFINE_TPLS(
  GTEST           "cmake/tpls/"  PT
  EIGEN           "cmake/tpls/"  PT
  PYBIND11	  "cmake/tpls/"  ST
  MKL             "cmake/tpls/"  EX
  quadmath        "cmake/tpls/"  EX
  BinUtils	  "TriBITS/tribits/common_tpls/"   ST
  BLAS            "TriBITS/tribits/common_tpls/"   PT
  LAPACK          "TriBITS/tribits/common_tpls/"   PT
  MPI             "cmake/tpls/"  ST
  TRILINOS        "cmake/tpls/"  ST
  BLAZE           "cmake/tpls/"  ST
  ARMADILLO       "cmake/tpls/"  ST
  MKL 		  "cmake/tpls/"  EX
  # we need to strip kokkos from trilinos but for now leave it
  # later on we would like to possibly build rompp wihtout
  # trilinos if needed but kokkos. So we need to have kokkos libs
  # separate from trilinos libraries
  #KOKKOS         "cmake/tpls/"  PT
  )

# PT means required
# ST optional
#MPI            "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/FindTPLMPI.cmake"   PT
