
add_serial_exe_and_test(integrate_to_time_with_failure
  ode integrate_to_target_time_with_failure.cc "PASSED")

add_serial_exe_and_test(arbitrary_stepper_1_aux_state_with_failure
  ode arbitrary_stepper_1_aux_state_with_failure.cc "PASSED")
add_serial_exe_and_test(arbitrary_stepper_2_aux_state_with_failure
  ode arbitrary_stepper_2_aux_state_with_failure.cc "PASSED")
add_serial_exe_and_test(arbitrary_stepper_3_aux_state_with_failure
  ode arbitrary_stepper_3_aux_state_with_failure.cc "PASSED")

add_serial_exe_and_test(bdf1_stepper_with_failure
  ode bdf1_stepper_with_failure.cc "PASSED")
# add_serial_exe_and_test(bdf2_stepper_with_failure
#   ode bdf2_stepper_with_failure.cc "PASSED")
