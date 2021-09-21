
set(TESTPREPENDEX ${TESTING_LEVEL}_ode_steppers_explicit)
set(TESTPREPENDIM ${TESTING_LEVEL}_ode_steppers_implicit)

add_serial_utest(${TESTPREPENDEX}_concepts ${CMAKE_CURRENT_SOURCE_DIR}/ode_concepts_explicit.cc)
add_serial_utest(${TESTPREPENDIM}_concepts ${CMAKE_CURRENT_SOURCE_DIR}/ode_concepts_implicit.cc)

if(PRESSIO_ENABLE_TPL_KOKKOS)
 	set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_stencil_data_kokkos.cc)
	add_serial_utest_kokkos(${TESTPREPENDIM}_stencil_data_kokkos ${SRC})
endif()

if(PRESSIO_ENABLE_TPL_EIGEN)
	set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_explicit_steppers_numerics_eigen.cc)
	add_serial_utest(${TESTPREPENDEX}_numerics_eigen ${SRC})

	set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_stencil_data_eigen.cc)
	add_serial_utest(${TESTPREPENDIM}_stencil_data_eigen ${SRC})

	set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_bdf1_numerics_eigen.cc)
	add_serial_utest(${TESTPREPENDIM}_bdf1_numerics_eigen ${SRC})

	set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_bdf2_numerics_eigen.cc)
	add_serial_utest(${TESTPREPENDIM}_bdf2_numerics_eigen ${SRC})

	set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_crank_nicolson_numerics_eigen.cc)
	add_serial_utest(${TESTPREPENDIM}_crank_nicolson_numerics_eigen ${SRC})

	set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_ensure_app_called_with_correct_time.cc)
	add_serial_utest(${TESTPREPENDIM}_ensure_app_called_with_correct_time_eigen ${SRC})
endif()

set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_explicit_steppers_numerics_custom_types.cc)
add_serial_utest(${TESTPREPENDEX}_numerics_custom_types ${SRC})

set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_steppers_numerics_custom_types.cc)
add_serial_utest(${TESTPREPENDIM}_numerics_custom_types ${SRC})

add_subdirectory(with_time_step_recovery)