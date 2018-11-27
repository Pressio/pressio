
add_serial_unit_test(config_defs_structs)
add_serial_unit_test(meta_basic)
add_serial_unit_test(meta_detect_operators)
add_serial_unit_test(meta_detect_typedefs)

add_serial_unit_test(vector_serial_eigen_dynamic_class)
add_serial_unit_test(vector_serial_eigen_dynamic_expr_templates)
add_serial_unit_test(matrix_vector_product_eigen)
add_serial_unit_test(matrix_matrix_product_eigen)
add_serial_unit_test(multi_vector_serial_eigen_dynamic_class)
add_serial_unit_test(multi_vector_eigen_product_eigen_vector)
add_serial_unit_test(sparse_matrix_product_multi_vector_eigen)
add_serial_unit_test(multi_vector_eigen_dot_eigen_vector)

IF(HAVE_BLAZE)
	add_serial_unit_test(vector_serial_blaze_class)
ENDIF()

IF(HAVE_ARMADILLO)
	add_serial_unit_test(vector_serial_armadillo_class)
ENDIF()


# SET(TESTFILES "")
# APPEND_SET(TESTFILES utest_vector_meta.cc)
# APPEND_SET(TESTFILES utest_vector_serial_eigen_traits.cc)
# APPEND_SET(TESTFILES utest_vector_serial_stdlib_traits.cc)
# APPEND_SET(TESTFILES utest_vector_serial_eigen_class.cc)
# APPEND_SET(TESTFILES utest_matrix_meta.cc)
# APPEND_SET(TESTFILES utest_matrix_dense_serial_eigen_traits.cc)
# APPEND_SET(TESTFILES utest_matrix_dense_serial_stdlib_traits.cc)
# APPEND_SET(TESTFILES utest_matrix_sparse_serial_eigen_traits.cc)
# APPEND_SET(TESTFILES utest_matrix_sparse_serial_eigen_class.cc)
# APPEND_SET(TESTFILES utest_matrix_dense_serial_eigen_class.cc)

# # append here tests to be grouped into a single executable
# # APPEND_SET(...)
# # ...
# # we create a single executable to run all the tests above
# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   unitTest_serial_grouped
#   SOURCES ${TESTFILES} gTestMain.cc
#   COMM serial mpi
#   NUM_MPI_PROCS 1
#   PASS_REGULAR_EXPRESSION "PASSED"
#   FAIL_REGULAR_EXPRESSION "FAILED"
# )
