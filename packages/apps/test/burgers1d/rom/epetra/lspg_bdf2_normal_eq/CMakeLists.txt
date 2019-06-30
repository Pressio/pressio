
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  burgers1d_rom_epetra_lspg_bdf2_gn_normal_eqns
  SOURCES main.cc
  NUM_MPI_PROCS 2
  PASS_REGULAR_EXPRESSION "PASSED"
  )

configure_file(../../basis/svd_basis_ncell20_t010_dt001_implicit_bdf2.txt basis.txt COPYONLY)
