TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_REQUIRED_PACKAGES mpl utils
  # LIB_OPTIONAL_PACKAGES
  #
  LIB_REQUIRED_TPLS EIGEN
  LIB_OPTIONAL_TPLS GTEST BLAS LAPACK MPI KOKKOS TRILINOS PYBIND11 BLAZE ARMADILLO BinUtils MKL
  #REGRESSION_EMAIL_LIST simplecxx-regressions@someurl.none
  )
