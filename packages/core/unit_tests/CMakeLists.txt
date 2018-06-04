

function(add_unit_test test_name ranks)
  TRIBITS_ADD_EXECUTABLE_AND_TEST(
    ${test_name} 
    SOURCES ${test_name}.cc gTestMain.cc
    COMM serial mpi
    NUM_MPI_PROCS ${ranks}
    STANDARD_PASS_OUTPUT
    #FAIL_REGULAR_EXPRESSION "  FAILED  "
  )
endfunction()


add_unit_test(Test1 1)
add_unit_test(Test2 2)


# function(add_unit_test test_name)
#   set(frontendtestfiles "${frontendtestfiles};${test_name}.cc" PARENT_SCOPE)
#   add_executable(${test_name} ${test_name}.cc gtest_main.cc)
#   # GMock links to GTest, so all we need is GMock
#   target_link_libraries(${test_name} ${DARMA_BACKEND_LIBNAME} ${GMOCK_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
#   add_test(${test_name} ${CMAKE_CURRENT_BINARY_DIR}/${test_name} "")
#   set_tests_properties(${test_name} PROPERTIES TIMEOUT 30
#     FAIL_REGULAR_EXPRESSION "FAILED;should be deleted but never is"
#     PASS_REGULAR_EXPRESSION "PASSED")
# endfunction()
