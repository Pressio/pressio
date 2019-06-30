
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  apps_unsteady_linear_adv_diff_1d_lspg_euler_gn_normal_eqns
  SOURCES main.cc
  NUM_MPI_PROCS 2
  PASS_REGULAR_EXPRESSION "PASSED"
  )

configure_file(../../basis/basis_row_39_cols_15.txt basis.txt COPYONLY)
