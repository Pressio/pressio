
set(ROOTNAME ${TESTING_LEVEL}_ode_implicit)

if(PRESSIO_ENABLE_TPL_EIGEN)

add_serial_utest(
  ${TESTING_LEVEL}_ode_implicit_arbitrary_stepper_run_eigen
  ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_arbitrary_stepper_correctness.cc)

add_serial_exe_and_test(
  burgers1d_bdf1_eigen ${ROOTNAME} burgers1d_bdf1_eigen.cc "PASSED")

add_serial_exe_and_test(
  burgers1d_bdf1_observe_eigen ${ROOTNAME} burgers1d_bdf1_observe_eigen.cc "PASSED")

add_serial_exe_and_test(
  burgers1d_bdf2_eigen ${ROOTNAME} burgers1d_bdf2_eigen.cc "PASSED")

add_serial_exe_and_test(
  burgers1d_bdf2_observe_eigen ${ROOTNAME} burgers1d_bdf2_observe_eigen.cc "PASSED")

add_serial_exe_and_test(
  burgers1d_cn_eigen ${ROOTNAME} burgers1d_cn_eigen.cc "PASSED")

endif()
