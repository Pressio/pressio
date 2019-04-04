
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  unsteady_linear_adv_diff_1d_fom_epetra_explicit_euler
  SOURCES main.cc
  NUM_MPI_PROCS 2
)
