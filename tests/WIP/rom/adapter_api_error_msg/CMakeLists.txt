
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

macro(add_serial_test_to_check_build_failure
    TESTNAME TESTSRCS GREPSTRING1 GREPSTRING2 PRIVATEDEF1 PRIVATEDEF2 PRIVATEDEF3)
  # set name of the executable
  set(testNameFinal rom_${TESTNAME})
  add_executable(${testNameFinal} ${TESTSRCS})

  # Avoid building these targets normally
  set_target_properties(${testNameFinal} PROPERTIES
    EXCLUDE_FROM_ALL TRUE
    EXCLUDE_FROM_DEFAULT_BUILD TRUE)

  target_compile_definitions(${testNameFinal}
    PRIVATE ${PRIVATEDEF1} ${PRIVATEDEF2} ${PRIVATEDEF3})

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
    -DstringToGrepOne=${GREPSTRING1}
    -DstringToGrepTwo=${GREPSTRING2}
    -P ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  # Expect these tests to fail (i.e. cmake --build should return a non-zero value)
  #set_tests_properties(${testNameFinal} PROPERTIES WILL_FAIL TRUE)

  # taken from: https://stackoverflow.com/questions/30155619/expected-build-failure-tests-in-cmake
endmacro()
#=====================================================================

set(dt_src	discrete_time_api.cc)
set(u_gale_src	galerkin_continuous_time_api.cc)
set(u_lspg_src	lspg_continuous_time_api.cc)
set(s_lspg_src  lspg_steady_api.cc)

set(LSPG_st "rom_compose_steady_lspg_impl")
set(LSPG_ust "rom_compose_unsteady_lspg_impl")
set(GALE "rom_compose_galerkin_impl")

set(MISS_VEL_C "is without (or has a wrong) create velocity method")
set(MISS_VEL   "is without (or has a wrong) velocity method")

set(MISS_APP_J_C  "is without (or has a wrong) create apply jacobian result method")
set(MISS_APP_J    "is without (or has a wrong) apply jacobian method")

set(MISS_RES_C "is without (or has a wrong) create residual method")
set(MISS_RES   "is without (or has a wrong) residual method")

set(MISS_DT_R_C "is without (or has a wrong) create discrete time residual method")
set(MISS_DT_R   "is without (or has a wrong) discrete time residual method")

set(MISS_APP_DT_J_C "is without (or has a wrong) create apply discrete time jacobian result method")
set(MISS_APP_DT_J   "is without (or has a wrong) apply discrete time jacobian method")

###############################
### DEFAULT LSPG STEADY ###
###############################
add_serial_test_to_check_build_failure(
  adapter_error_default_lspg_steady_api_case1
  ${s_lspg_src} ${MISS_RES} ${LSPG_st} NO_RES DEFAULT DO_LSPG)

add_serial_test_to_check_build_failure(
  adapter_error_default_lspg_steady_api_case2
  ${s_lspg_src} ${MISS_APP_J} ${LSPG_st} NO_APP_J DEFAULT DO_LSPG)

add_serial_test_to_check_build_failure(
  adapter_error_default_lspg_steady_api_case3
  ${s_lspg_src} ${MISS_RES_C} ${LSPG_st} NO_RES_C DEFAULT DO_LSPG)

add_serial_test_to_check_build_failure(
  adapter_error_default_lspg_steady_api_case4
  ${s_lspg_src} ${MISS_APP_J_C} ${LSPG_st} NO_APP_J_C DEFAULT DO_LSPG)

##################################
### PRECONDITIONED LSPG STEADY ###
##################################
add_serial_test_to_check_build_failure(
  adapter_error_precond_lspg_steady_api_case1
  ${s_lspg_src} ${MISS_RES} ${LSPG_st} NO_RES PRECOND DO_LSPG)

add_serial_test_to_check_build_failure(
  adapter_error_precond_lspg_steady_api_case2
  ${s_lspg_src} ${MISS_APP_J} ${LSPG_st} NO_APP_J PRECOND DO_LSPG)

add_serial_test_to_check_build_failure(
  adapter_error_precond_lspg_steady_api_case3
  ${s_lspg_src} ${MISS_RES_C} ${LSPG_st} NO_RES_C PRECOND DO_LSPG)

add_serial_test_to_check_build_failure(
  adapter_error_precond_lspg_steady_api_case4
  ${s_lspg_src} ${MISS_APP_J_C} ${LSPG_st} NO_APP_J_C PRECOND DO_LSPG)

##################################
### DEFAULT LSPG CONT TIME API ###
##################################
add_serial_test_to_check_build_failure(
  adapter_error_default_lspg_continuous_time_api_case1
  ${u_lspg_src} ${MISS_VEL} ${LSPG_ust} NO_VEL DEFAULT DO_LSPG)

add_serial_test_to_check_build_failure(
  adapter_error_default_lspg_continuous_time_api_case2
  ${u_lspg_src} ${MISS_APP_J} ${LSPG_ust} NO_APP_J DEFAULT DO_LSPG)

add_serial_test_to_check_build_failure(
  adapter_error_default_lspg_continuous_time_api_case3
  ${u_lspg_src} ${MISS_VEL_C} ${LSPG_ust} NO_VEL_C DEFAULT DO_LSPG)

add_serial_test_to_check_build_failure(
  adapter_error_default_lspg_continuous_time_api_case4
  ${u_lspg_src} ${MISS_APP_J_C} ${LSPG_ust} NO_APP_J_C DEFAULT DO_LSPG)


##########################################
### DEFAULT LSPG DISCRETE TIME API ###
##########################################
add_serial_test_to_check_build_failure(
  adapter_error_default_lspg_discrete_time_api_case1
  ${dt_src} ${MISS_DT_R} ${LSPG_ust} NO_RES DEFAULT DO_LSPG)

add_serial_test_to_check_build_failure(
  adapter_error_default_lspg_discrete_time_api_case2
  ${dt_src} ${MISS_APP_DT_J} ${LSPG_ust} NO_APP_J DEFAULT DO_LSPG)

add_serial_test_to_check_build_failure(
  adapter_error_default_lspg_discrete_time_api_case3
  ${dt_src} ${MISS_DT_R_C} ${LSPG_ust} NO_RES_C DEFAULT DO_LSPG)

add_serial_test_to_check_build_failure(
  adapter_error_default_lspg_discrete_time_api_case4
  ${dt_src} ${MISS_APP_DT_J_C} ${LSPG_ust} NO_APP_J_C DEFAULT DO_LSPG)



# ##################################
# ### DEFAULT GALERKIN CONT TIME ###
# ##################################
# add_serial_test_to_check_build_failure(
#   adapter_error_default_galerkin_continuous_time_api_case1
#   ${u_gale_src} ${MISS_VEL} ${GALE} NO_VEL DEFAULT DO_GALE)

# add_serial_test_to_check_build_failure(
#   adapter_error_default_galerkin_continuous_time_api_case2
#   ${u_gale_src} ${MISS_VEL_C} ${GALE} NO_VEL_C DEFAULT DO_GALE)

# ##########################################
# ### DEFAULT GALERKIN DISCRETE TIME API ###
# ##########################################
# add_serial_test_to_check_build_failure(
#   adapter_error_default_galerkin_discrete_time_api_case1
#   ${dt_src} ${MISS_DT_R} ${GALE} NO_RES DEFAULT DO_GALE)

# add_serial_test_to_check_build_failure(
#   adapter_error_default_galerkin_discrete_time_api_case2
#   ${dt_src} ${MISS_APP_DT_J} ${GALE} NO_APP_J DEFAULT DO_GALE)

# add_serial_test_to_check_build_failure(
#   adapter_error_default_galerkin_discrete_time_api_case3
#   ${dt_src} ${MISS_DT_R_C} ${GALE} NO_RES_C DEFAULT DO_GALE)

# add_serial_test_to_check_build_failure(
#   adapter_error_default_galerkin_discrete_time_api_case4
#   ${dt_src} ${MISS_APP_DT_J_C} ${GALE} NO_APP_J_C DEFAULT DO_GALE)

# add_serial_test_to_check_build_failure(
#   adapter_error_default_galerkin_discrete_time_api_case5
#   ${dt_src} "is without (or has a wrong) dense matrix typedef" ${GALE} NO_MAT_ALIAS DEFAULT DO_GALE)

#add_serial_exe_and_test(discrete_api_error_msg rom discrete_api_error_msg.cc "PASSED")
#add_serial_exe_and_test(continuous_time_api_error_msg rom continuous_time_api_error_msg.cc "PASSED")
