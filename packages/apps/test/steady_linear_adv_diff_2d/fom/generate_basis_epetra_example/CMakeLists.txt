
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  steadyAdvDiff2d_fom_gen_basis_epetra_example
  SOURCES main.cc
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  )

configure_file(../../rom/basis/basis_row189_col5.txt gold_basis.txt COPYONLY)
