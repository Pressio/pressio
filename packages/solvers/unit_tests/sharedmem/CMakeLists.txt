

SET(TESTFILES "")

APPEND_SET(TESTFILES test_solvers_linear_dense_eigen.cc)
APPEND_SET(TESTFILES test_solvers_linear_eigen.cc)
APPEND_SET(TESTFILES test_system_traits.cc)
APPEND_SET(TESTFILES test_solvers_factory_nonlinear.cc)
#APPEND_SET(TESTFILES test_solvers_nonlinear_policies.cc)
APPEND_SET(TESTFILES test_solvers_nonlinear_eigen.cc)
APPEND_SET(TESTFILES test_solvers_leastsquare_linear_eigen.cc)
APPEND_SET(TESTFILES test_solvers_nonlinear_leastsquare_eigen.cc)

# append here tests to be grouped into a single executable
# APPEND_SET(...)
# ...

# we create a single executable to run all the tests above
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  solvers_unitTest_serial_grouped
  SOURCES ${TESTFILES} gTestMain.cc
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  FAIL_REGULAR_EXPRESSION "FAILED"
)


TRIBITS_ADD_EXECUTABLE_AND_TEST(
  gauss_newton_with_qr_inner_solver
  SOURCES gauss_newton_with_qr_inner_solver gTestMain.cc
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  FAIL_REGULAR_EXPRESSION "FAILED"
)
