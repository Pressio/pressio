cmake_minimum_required(VERSION 3.11)

TRIBITS_ADD_EXECUTABLE_AND_TEST(
  laplace1d_epetra_fom
  SOURCES ../laplace1dEpetra.cc main.cc
  NUM_MPI_PROCS 3
  #STANDARD_PASS_OUTPUT
  #PASS_REGULAR_EXPRESSION "Hello World"
  )
