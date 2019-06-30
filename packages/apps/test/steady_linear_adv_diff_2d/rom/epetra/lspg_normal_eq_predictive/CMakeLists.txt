
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  steady_lin_adv_diff2d_rom_epetra_lspg_gn_normal_eqns_predictive
  SOURCES main.cc
  NUM_MPI_PROCS 3
  PASS_REGULAR_EXPRESSION "PASSED"
  )

configure_file(../../basis/basis_row189_col5.txt basis.txt COPYONLY)
