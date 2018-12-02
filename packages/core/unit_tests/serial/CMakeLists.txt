

ADD_SUBDIRECTORY(meta)
ADD_SUBDIRECTORY(vector)
ADD_SUBDIRECTORY(multi_vector)
# ADD_SUBDIRECTORY(matrix)
ADD_SUBDIRECTORY(dot)
ADD_SUBDIRECTORY(product)
ADD_SUBDIRECTORY(miscellaneous)


# SET(TESTFILES "")
# APPEND_SET(TESTFILES utest_vector_meta.cc)
# APPEND_SET(TESTFILES utest_vector_serial_eigen_traits.cc)

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
