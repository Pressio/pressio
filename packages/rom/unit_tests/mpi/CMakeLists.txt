

# this function creates and adds an individual test executable
# this is needed for cases where we need to specify the number
# of MPI processes, for example to do testing for distributed objects
function(add_single_unit_test test_name ranks)  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    rom_${test_name}
    SOURCES ${test_name}.cc gTestMain.cc
    COMM serial mpi
    NUM_MPI_PROCS ${ranks}
    PASS_REGULAR_EXPRESSION "PASSED"
    FAIL_REGULAR_EXPRESSION "FAILED"
   )
endfunction()

# append tests need to have a unique executable each
# add_single_unit_test(utest_matrix_pseudoInverse_distributed_epetra 3)

