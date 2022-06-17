
set(TESTPREPEND ${TESTING_LEVEL}_ode_steppers_explicit)

add_serial_utest(${TESTPREPEND}_concepts ${CMAKE_CURRENT_SOURCE_DIR}/ode_concepts_explicit.cc)

if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_explicit_steppers_numerics_eigen.cc)
  add_serial_utest(${TESTPREPEND}_numerics_eigen ${SRC})

  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_euler_mass_matrix_numerics_eigen.cc)
  add_serial_utest(${TESTPREPEND}_euler_mass_matrix_numerics_eigen ${SRC})
endif()

set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_explicit_steppers_numerics_custom_types.cc)
add_serial_utest(${TESTPREPEND}_numerics_custom_types ${SRC})
