
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  nonlin_unsteady_adv_diff_reac2d_eigen_gen_basis_implicit_bdf2
  SOURCES main.cc
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  )

configure_file(../../rom/basis/basis_nx11_ny21_col11_bdf2.txt gold_basis.txt COPYONLY)
