
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  burgers1d_epetra_lspg_euler_gn_nor_eqns_reduced_no_masked
  SOURCES ../burgers1dEpetra.cc ../burgers1dEpetra_masked.cc main_gn_nor_eqns.cc
  NUM_MPI_PROCS 2
  #STANDARD_PASS_OUTPUT
  #PASS_REGULAR_EXPRESSION "Hello World"
  )

configure_file(../U20_ncell50_t35_dt001.txt U20_ncell50_t35_dt001.txt COPYONLY)
