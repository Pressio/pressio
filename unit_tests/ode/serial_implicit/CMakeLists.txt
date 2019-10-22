
add_serial_utest(ode_model_admissilble_implicit_stepper
  ${CMAKE_CURRENT_SOURCE_DIR}/ode_model_admissilble_implicit_stepper.cc)

add_serial_utest(ode_implicit_euler ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_euler.cc)
add_serial_utest(ode_implicit_bdf2 ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_bdf2.cc)

add_serial_utest(ode_implicit_arbitrary_stepper
  ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_arbitrary_stepper.cc)
