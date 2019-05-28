
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  advDiffReacFlame2d_fom_smeigen_explicit_rk4
  SOURCES main.cc
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  )

configure_file(../plot_help.py plot_help.py COPYONLY)
configure_file(./mesh.dat mesh.dat COPYONLY)
configure_file(./gappy_to_full_gid_map.dat gappy_to_full_gid_map.dat COPYONLY)
