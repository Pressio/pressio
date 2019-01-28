
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  burgers1d_epetra_lspg_bdf2_gn_nor_eqns
  SOURCES ../burgers1d_epetra_lspg_euler/burgers1dEpetra.cc main_gn_nor_eqns.cc
  NUM_MPI_PROCS 2
  #STANDARD_PASS_OUTPUT
  #PASS_REGULAR_EXPRESSION "Hello World"
  )

# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   burgers1d_epetra_lspg_euler_gn_qr
#   SOURCES main_gn_qr.cc
#   NUM_MPI_PROCS 2
#   #STANDARD_PASS_OUTPUT
#   #PASS_REGULAR_EXPRESSION "Hello World"
#   )

configure_file(U20_ncell50_t35_dt001.txt U20_ncell50_t35_dt001.txt COPYONLY)
