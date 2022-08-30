
set(TESTPREPEND ${TESTING_LEVEL}_ode_steppers_implicit)

if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_ensure_app_called_with_correct_time.cc)
  add_serial_utest(${TESTPREPEND}_ensure_app_called_with_correct_time_eigen ${SRC})
endif()

set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_steppers_numerics_custom_types.cc)
add_serial_utest(${TESTPREPEND}_numerics_custom_types ${SRC})

# add_subdirectory(with_time_step_recovery)
