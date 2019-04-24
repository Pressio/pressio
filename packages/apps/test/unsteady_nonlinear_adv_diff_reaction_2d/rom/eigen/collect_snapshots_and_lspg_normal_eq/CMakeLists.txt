
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  nonlin_unsteady_adv_diff_reac2d_eigen_do_snaps_and_lspg_normal_eqns
  SOURCES main.cc
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  )

configure_file(../../basis/basis_nx11_ny21_col11.txt gold_basis.txt COPYONLY)
