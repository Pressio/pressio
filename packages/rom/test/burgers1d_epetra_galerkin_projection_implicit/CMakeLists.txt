TRIBITS_ADD_EXECUTABLE_AND_TEST(
  rom_burgers1d_epetra_galpro_implicit
  SOURCES main.cc
  NUM_MPI_PROCS 3
  #STANDARD_PASS_OUTPUT
  )
