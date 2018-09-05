


# # this function creates and adds an individual test executable
# function(add_single_unit_test test_name)  
#   TRIBITS_ADD_EXECUTABLE_AND_TEST(
#     ${test_name}
#     SOURCES ${test_name}.cc gTestMain.cc
#     COMM serial mpi
#     PASS_REGULAR_EXPRESSION "PASSED"
#     FAIL_REGULAR_EXPRESSION "FAILED"
#   )
# endfunction()

# add_single_unit_test(utest_vector_distributed_epetra_traits 1)



SET(TESTFILES "")

APPEND_SET(TESTFILES utest_vector_sharedmem_kokkos_class.cc)

# append here tests to be grouped into a single executable
# APPEND_SET(...)
# ...
# we create a single executable to run all the tests above
TRIBITS_ADD_EXECUTABLE_AND_TEST(
  unitTest_kokkos_grouped
  SOURCES ${TESTFILES} gTestMain.cc
  COMM serial mpi
  NUM_MPI_PROCS 1
  PASS_REGULAR_EXPRESSION "PASSED"
  FAIL_REGULAR_EXPRESSION "FAILED"
)
