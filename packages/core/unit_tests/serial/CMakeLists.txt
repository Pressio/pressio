
add_unit_test(meta_basic ${GTESTMAINSDIR}/gTestMain_serial)
add_unit_test(meta_detect_operators ${GTESTMAINSDIR}/gTestMain_serial)
add_unit_test(meta_detect_typedefs ${GTESTMAINSDIR}/gTestMain_serial)


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
# APPEND_SET(TESTFILES utest_matrix_vector_product.cc)

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
