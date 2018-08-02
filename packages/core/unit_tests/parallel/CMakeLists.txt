

# this function creates and adds an individual test executable
# this is needed for cases where we need to specify the number
# of MPI processes, for example to do testing for distributed objects
function(add_single_unit_test test_name ranks)  
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    core_${test_name}
    SOURCES ${test_name}.cc gTestMain.cc
    COMM serial mpi
    NUM_MPI_PROCS ${ranks}
    PASS_REGULAR_EXPRESSION "PASSED"
    FAIL_REGULAR_EXPRESSION "FAILED"
  )
endfunction()

# append tests need to have a unique executable each
add_single_unit_test(test_core_vector_distributed_epetra_traits 1)
add_single_unit_test(test_core_vector_distributed_epetra_class 3)
add_single_unit_test(test_core_matrix_vector_product 1)
add_single_unit_test(test_core_vector_dot_product 1)
add_single_unit_test(test_core_matrix_sparse_distributed_epetra_class 3)
add_single_unit_test(test_core_matrix_matrix_product_sparse_distributed_epetra 3)
add_single_unit_test(test_core_multi_vector_distributed_epetra_traits 1)
add_single_unit_test(test_core_multi_vector_distributed_epetra_class 3)
add_single_unit_test(test_core_matrix_dense_distributed_epetra_traits 3)
add_single_unit_test(test_core_matrix_dense_distributed_epetra_class 3)

add_single_unit_test(test_core_matrix_matrix_product_dense_distributed_epetra 3)
add_single_unit_test(test_core_matrix_matrix_product_dense_sparse_distributed_epetra 3)
