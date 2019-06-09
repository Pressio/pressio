
TRIBITS_ADD_EXECUTABLE(
  ks1d_fom_eigen_implicit_bdf2
  SOURCES main.cc
  )

ADD_SUBDIRECTORIES(transient predictive reproductive)