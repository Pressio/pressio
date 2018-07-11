

SET(TESTFILES "")

#APPEND_SET(TESTFILES test_solvers_linear_eigen.cc)
#APPEND_SET(TESTFILES test_solvers_nonlinear_eigen.cc)
APPEND_SET(TESTFILES test_solvers_linear_iterative_factory.cc)

# append here tests to be grouped into a single executable
# APPEND_SET(...)
# ...


# we create a single executable to run all the tests above
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  solvers_unitTest_serial_grouped
  SOURCES ${TESTFILES} gTestMain.cc
  COMM serial mpi
  NUM_MPI_PROCS 1
  FAIL_REGULAR_EXPRESSION "FAILED"
)
