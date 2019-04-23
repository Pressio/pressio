
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  nonlin_unsteady_adv_diff_reac2d_epetra_do_snaps_and_lspg_normal_eqns
  SOURCES main.cc
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  )
