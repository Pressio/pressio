
add_serial_utest(
	${TESTING_LEVEL}_ode_concepts_common
	${CMAKE_CURRENT_SOURCE_DIR}/ode_concepts_common.cc)

add_serial_utest(
	${TESTING_LEVEL}_ode_concepts_explicit
	${CMAKE_CURRENT_SOURCE_DIR}/ode_concepts_explicit.cc)

add_serial_utest(
	${TESTING_LEVEL}_ode_concepts_implicit
	${CMAKE_CURRENT_SOURCE_DIR}/ode_concepts_implicit.cc)

add_serial_utest(
	${TESTING_LEVEL}_ode_advance_functions
	${CMAKE_CURRENT_SOURCE_DIR}/ode_advance_functions.cc)

if(PRESSIO_ENABLE_TPL_EIGEN)
	set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_explicit_numerics_eigen.cc)
	add_serial_utest(${TESTING_LEVEL}_ode_explicit_steppers_numerics_eigen ${SRC})

	set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_numerics_eigen.cc)
	add_serial_utest(${TESTING_LEVEL}_ode_implicit_steppers_numerics_eigen ${SRC})

	set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_crank_nic_numerics_eigen.cc)
	add_serial_utest(${TESTING_LEVEL}_ode_implicit_crank_nicolson_stepper_numerics_eigen ${SRC})

	set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_implicit_app_called_with_correct_time.cc)
	add_serial_utest(${TESTING_LEVEL}_ode_implicit_app_called_with_correct_time_eigen ${SRC})
endif()

set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/ode_explicit_numerics_custom_types.cc)
add_serial_utest(${TESTING_LEVEL}_ode_explicit_steppers_numerics_custom_types ${SRC})


add_subdirectory(advance_with_time_step_recovery)