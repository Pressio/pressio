
TRIBITS_ADD_EXECUTABLE(
  steadyAdvDiff2d_fom_epetra
  SOURCES main.cc
  )

TRIBITS_ADD_TEST(
  steadyAdvDiff2d_fom_epetra  NUM_MPI_PROCS 1
  )

TRIBITS_ADD_TEST(
  steadyAdvDiff2d_fom_epetra  NUM_MPI_PROCS 2
  )

TRIBITS_ADD_TEST(
  steadyAdvDiff2d_fom_epetra  NUM_MPI_PROCS 3
  )
