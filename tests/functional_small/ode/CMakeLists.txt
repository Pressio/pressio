
add_serial_utest(
	${TESTING_LEVEL}_ode_concepts_advancers
	${CMAKE_CURRENT_SOURCE_DIR}/ode_concepts_advancers.cc)

add_serial_utest(
	${TESTING_LEVEL}_ode_concepts_explicit
	${CMAKE_CURRENT_SOURCE_DIR}/ode_concepts_explicit.cc)

add_serial_utest(
	${TESTING_LEVEL}_ode_concepts_implicit
	${CMAKE_CURRENT_SOURCE_DIR}/ode_concepts_implicit.cc)

add_serial_utest(
	${TESTING_LEVEL}_ode_advance_functions
	${CMAKE_CURRENT_SOURCE_DIR}/ode_advance_functions.cc)

if(PRESSIO_ENABLE_TPL_KOKKOS)
 	set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_stencil_data_kokkos.cc)
	add_serial_utest_kokkos(${TESTING_LEVEL}_ode_implicit_stencil_data_kokkos ${SRC})
endif()

if(PRESSIO_ENABLE_TPL_EIGEN)
	set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_explicit_steppers_numerics_eigen.cc)
	add_serial_utest(${TESTING_LEVEL}_ode_explicit_steppers_numerics_eigen ${SRC})

	set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_stencil_data_eigen.cc)
	add_serial_utest(${TESTING_LEVEL}_ode_implicit_stencil_data_eigen ${SRC})

	set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_bdf1_numerics_eigen.cc)
	add_serial_utest(${TESTING_LEVEL}_ode_implicit_bdf1_numerics_eigen ${SRC})

	set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_bdf2_numerics_eigen.cc)
	add_serial_utest(${TESTING_LEVEL}_ode_implicit_bdf2_numerics_eigen ${SRC})

	set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_crank_nicolson_numerics_eigen.cc)
	add_serial_utest(${TESTING_LEVEL}_ode_implicit_crank_nicolson_stepper_numerics_eigen ${SRC})

	set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_ensure_app_called_with_correct_time.cc)
	add_serial_utest(${TESTING_LEVEL}_ode_implicit_ensure_app_called_with_correct_time_eigen ${SRC})
endif()

set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_explicit_steppers_numerics_custom_types.cc)
add_serial_utest(${TESTING_LEVEL}_ode_explicit_steppers_numerics_custom_types ${SRC})


add_subdirectory(advance_with_time_step_recovery)