TRIBITS_REPOSITORY_DEFINE_TPLS(
  EIGEN           "cmake/tpls/"  PT
  GTEST           "cmake/tpls/"  ST
  PYBIND11	  "cmake/tpls/"  ST
  BinUtils	  "TriBITS/tribits/common_tpls/"   ST
  BLAS            "TriBITS/tribits/common_tpls/"   ST
  LAPACK          "TriBITS/tribits/common_tpls/"   ST
  MPI             "cmake/tpls/"  ST
  KOKKOS	  "cmake/tpls/"  ST
  TRILINOS        "cmake/tpls/"  ST
  BLAZE           "cmake/tpls/"  ST
  ARMADILLO       "cmake/tpls/"  ST
  MKL             "cmake/tpls/"  EX
  quadmath        "cmake/tpls/"  EX
  )

# PT means required
# ST optional
#MPI            "${${PROJECT_NAME}_TRIBITS_DIR}/core/std_tpls/FindTPLMPI.cmake"   PT
