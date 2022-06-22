
set(TESTPREPEND ${TESTING_LEVEL}_ode_steppers_implicit)

if(PRESSIO_ENABLE_TPL_KOKKOS)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_stencil_data_kokkos.cc)
  add_serial_utest_kokkos(${TESTPREPEND}_stencil_data_kokkos ${SRC})
endif()

if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_stencil_data_eigen.cc)
  add_serial_utest(${TESTPREPEND}_stencil_data_eigen ${SRC})

  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_bdf1_numerics_eigen.cc)
  add_serial_utest(${TESTPREPEND}_bdf1_numerics_eigen ${SRC})

  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_bdf2_numerics_eigen.cc)
  add_serial_utest(${TESTPREPEND}_bdf2_numerics_eigen ${SRC})

  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_crank_nicolson_numerics_eigen.cc)
  add_serial_utest(${TESTPREPEND}_crank_nicolson_numerics_eigen ${SRC})

  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_ensure_app_called_with_correct_time.cc)
  add_serial_utest(${TESTPREPEND}_ensure_app_called_with_correct_time_eigen ${SRC})
endif()

set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_steppers_numerics_custom_types.cc)
add_serial_utest(${TESTPREPEND}_numerics_custom_types ${SRC})

# add_subdirectory(with_time_step_recovery)
