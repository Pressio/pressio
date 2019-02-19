
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  burgers1d_epetra_galerkin_explicit_euler
  SOURCES ../burgers1dEpetra.cc main.cc
  NUM_MPI_PROCS 2
  )

configure_file(U20_ncell50_t35_dt001.txt U20_ncell50_t35_dt001.txt COPYONLY)
