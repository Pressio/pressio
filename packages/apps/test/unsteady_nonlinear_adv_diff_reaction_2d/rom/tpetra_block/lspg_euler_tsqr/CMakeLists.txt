
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  nonlin_unsteady_adv_diff_reac2d_tpetra_block_lspg_euler_tsqr
  SOURCES main.cc
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  )

configure_file(../../basis/basis_nx11_ny21_col11.txt basis.txt COPYONLY)
