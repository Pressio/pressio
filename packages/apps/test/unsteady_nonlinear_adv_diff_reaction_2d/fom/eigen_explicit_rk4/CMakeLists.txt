
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  advDiffReac2d_fom_eigen_explicit_rk4
  SOURCES main.cc
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  )

configure_file(../plot_help.py plot_help.py COPYONLY)
