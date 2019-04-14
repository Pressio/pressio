
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  steady_adv_diff1d_rom_epetra_lspg_gn_qr
  SOURCES main.cc
  NUM_MPI_PROCS 2
  PASS_REGULAR_EXPRESSION "PASSED"
  )

configure_file(../../basis/basis_row199_col15.txt basis.txt COPYONLY)
