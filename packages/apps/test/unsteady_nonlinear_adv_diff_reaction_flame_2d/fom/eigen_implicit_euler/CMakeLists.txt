
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  advDiffReacFlame2d_fom_eigen_implicit_euler
  SOURCES main.cc
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  )

configure_file(../plot_help.py plot_help.py COPYONLY)
