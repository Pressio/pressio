
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  advDiffReacFlame2d_eigen_do_snaps_and_lspg_normal_eqns_mask
  SOURCES main.cc
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  )

configure_file(r_gids.dat r_gids.dat COPYONLY)
