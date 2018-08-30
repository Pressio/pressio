TRIBITS_ADD_EXECUTABLE_AND_TEST(
  epetra_wrapper
  SOURCES main.cc
  NUM_MPI_PROCS 2
  #STANDARD_PASS_OUTPUT
  )

configure_file(mat_svd.dat mat_svd.dat COPYONLY)
configure_file(mat_svd_notsq.dat mat_svd_notsq.dat COPYONLY)
