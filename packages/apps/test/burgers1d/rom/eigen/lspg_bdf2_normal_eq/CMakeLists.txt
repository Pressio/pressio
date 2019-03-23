
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  burgers1d_rom_eigen_lspg_bdf2_gn_normal_eqns
  SOURCES main.cc
  NUM_MPI_PROCS 1
  #STANDARD_PASS_OUTPUT
  #PASS_REGULAR_EXPRESSION "Hello World"
  )

configure_file(../../basis/svd_basis_ncell20_t010_dt001_implicit_bdf2.txt basis.txt COPYONLY)
