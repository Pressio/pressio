
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  burgers1d_epetra_galerkin_explicit_euler
  SOURCES main.cc
  NUM_MPI_PROCS 2
  PASS_REGULAR_EXPRESSION "PASSED"
  )

configure_file(U20_ncell50_t35_dt001.txt basis.txt COPYONLY)
