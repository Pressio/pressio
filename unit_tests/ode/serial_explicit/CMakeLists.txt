
add_serial_utest(ode_serial_explicit_meta
  ${CMAKE_CURRENT_SOURCE_DIR}/ode_user_defined_model_admissibility.cc)

set(SRC2
  ${CMAKE_CURRENT_SOURCE_DIR}/ode_explicit_euler_user_defined_ops.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ode_explicit_euler.cc)
add_serial_utest(ode_serial_explicit_euler "${SRC2}")

set(SRC3
  ${CMAKE_CURRENT_SOURCE_DIR}/ode_explicit_rk4_user_defined_ops.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/ode_explicit_rk4.cc)
add_serial_utest(ode_serial_explicit_rk4 "${SRC3}")
