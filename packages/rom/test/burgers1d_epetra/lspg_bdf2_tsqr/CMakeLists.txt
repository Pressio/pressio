
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  burgers1d_epetra_lspg_bdf2_tsqr
  SOURCES ../burgers1dEpetra.cc main.cc
  NUM_MPI_PROCS 2
  #STANDARD_PASS_OUTPUT
  #PASS_REGULAR_EXPRESSION "Hello World"
  )

configure_file(../U20_ncell50_t35_dt001.txt U20_ncell50_t35_dt001.txt COPYONLY)
