
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  steady_adv_diff1d_rom_epetra_lspg_gn_normal_eqns
  SOURCES main.cc
  NUM_MPI_PROCS 2
  #STANDARD_PASS_OUTPUT
  #PASS_REGULAR_EXPRESSION "Hello World"
  )

configure_file(../../basis/fake_basis_row20_col5.txt basis.txt COPYONLY)
