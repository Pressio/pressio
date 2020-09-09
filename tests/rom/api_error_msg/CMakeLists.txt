
# the following macro creates and adds an individual SERIAL test that
# we use for testing that the build of an exe fails.
# The way this works is that the test itself is excluded from
# all the list of tests built regularly and we make
# ctest itself try to build the exe.
# When we run ctest, if the build fails, it means the test succeeds.

# WHy? Suppose we want to test that when the user has a wrong
# adapter class, they cannot create a ROM problem for a particular reason
# which is triggered by pressio at compile time.
# We want to test that for specific scenarios one cannot even build the executable.

macro(add_serial_test_to_check_build_failure TESTNAME PKGNAME TESTSRCS PRIVATEDEFS GREPSTRING)
  # set name of the executable
  set(testNameFinal ${PKGNAME}_${TESTNAME})
  add_executable(${testNameFinal} ${TESTSRCS})

  # Avoid building these targets normally
  set_target_properties(${testNameFinal} PROPERTIES
    EXCLUDE_FROM_ALL TRUE
    EXCLUDE_FROM_DEFAULT_BUILD TRUE)

  target_compile_definitions(${testNameFinal} PRIVATE ${PRIVATEDEFS})

  # Add the tests.  These invoke "cmake --build ..." which is a
  # cross-platform way of building the given target.
  # add_test(
  #   NAME ${testNameFinal}
  #   COMMAND ${CMAKE_COMMAND} --build . --target ${testNameFinal} --config $<CONFIGURATION>
  #   WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  # here we use runtest.cmake which will
  # - attempt to buold the executable for the test
  # - dump the log to file
  # - check that the file contains specific string
  add_test(
    NAME ${testNameFinal}
    COMMAND ${CMAKE_COMMAND}
    -DtargetName=${testNameFinal}
    -DtargetLogFile=${CMAKE_CURRENT_BINARY_DIR}/log_${testNameFinal}.out
    -DstringToGrep=${GREPSTRING}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  # Expect these tests to fail (i.e. cmake --build should return a non-zero value)
  #set_tests_properties(${testNameFinal} PROPERTIES WILL_FAIL TRUE)

  # taken from: https://stackoverflow.com/questions/30155619/expected-build-failure-tests-in-cmake
endmacro()
#=====================================================================

#####################
### DISCRETE API ###
#####################
add_serial_test_to_check_build_failure(
  discrete_api_error_fail1 rom discrete_api_error_msg.cc
  FAILMODE1 "missing the discrete_time_residual_method")

add_serial_test_to_check_build_failure(
  discrete_api_error_fail2 rom discrete_api_error_msg.cc
  FAILMODE2 "missing the apply_discrete_time_jacobian_method")

add_serial_test_to_check_build_failure(
  discrete_api_error_fail3 rom discrete_api_error_msg.cc
  FAILMODE3 "missing the create_discrete_time_residual_method")

add_serial_test_to_check_build_failure(
  discrete_api_error_fail4 rom discrete_api_error_msg.cc
  FAILMODE4 "missing the create_apply_discrete_time_jacobian_result_method")

add_serial_test_to_check_build_failure(
  discrete_api_error_fail5 rom discrete_api_error_msg.cc
  FAILMODE5 "missing the dense_matrix typedef")

#####################
### CONT TIME API ###
#####################
add_serial_test_to_check_build_failure(
  continuous_time_api_error_fail1 rom continuous_time_api_error_msg.cc
  FAILMODE1 "missing the velocity_method")

add_serial_test_to_check_build_failure(
  continuous_time_api_error_fail2 rom continuous_time_api_error_msg.cc
  FAILMODE2 "missing the apply_jacobian_method")

add_serial_test_to_check_build_failure(
  continuous_time_api_error_fail3 rom continuous_time_api_error_msg.cc
  FAILMODE3 "missing the create_velocity_method")

add_serial_test_to_check_build_failure(
  continuous_time_api_error_fail4 rom continuous_time_api_error_msg.cc
  FAILMODE4 "missing the create_apply_jacobian_result_method")


#add_serial_exe_and_test(discrete_api_error_msg rom discrete_api_error_msg.cc "PASSED")
#add_serial_exe_and_test(continuous_time_api_error_msg rom continuous_time_api_error_msg.cc "PASSED")
