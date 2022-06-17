
set(ROOTNAME ${TESTING_LEVEL}_ode_steppers_implicit)

if(PRESSIO_ENABLE_TPL_EIGEN)
add_serial_exe_and_test(
  integrate_to_time_with_failure_arbitrary_stepper_1_aux_state
  ${ROOTNAME} arbitrary_stepper_1_aux_state.cc "PASSED")

add_serial_exe_and_test(
  integrate_to_time_with_failure_arbitrary_stepper_2_aux_state
  ${ROOTNAME} arbitrary_stepper_2_aux_state.cc "PASSED")

add_serial_exe_and_test( 
  integrate_to_time_with_failure_arbitrary_stepper_3_aux_state
  ${ROOTNAME} arbitrary_stepper_3_aux_state.cc "PASSED")

add_serial_exe_and_test(
  integrate_to_time_with_failure_bdf1_stepper
  ${ROOTNAME} bdf1_stepper.cc "PASSED")

endif()
