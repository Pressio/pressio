
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  steady_adv_diff2d_rom_epetra_lspg_gn_normal_eqns_predictive
  SOURCES main.cc
  NUM_MPI_PROCS 3
  #STANDARD_PASS_OUTPUT
  #PASS_REGULAR_EXPRESSION "Hello World"
  )

configure_file(../../basis/basis_row189_col5.txt basis.txt COPYONLY)
