TRIBITS_ADD_EXECUTABLE_AND_TEST(
  burgers1d_fom_eigen_implicit_euler_observer
  SOURCES main.cc
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  )
