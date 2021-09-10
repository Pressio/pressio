
set(ROOTNAME ${TESTING_LEVEL}_ode_advancers)

add_serial_utest(
	${ROOTNAME}_concepts
	${CMAKE_CURRENT_SOURCE_DIR}/ode_concepts_advancers.cc)

add_serial_utest(
	${ROOTNAME}_functions 
	${CMAKE_CURRENT_SOURCE_DIR}/ode_advance_functions.cc)

if(PRESSIO_ENABLE_TPL_EIGEN)
add_serial_exe_and_test(
	integrate_to_time_with_failure_mock_stepper
	${ROOTNAME} to_target_time_with_time_step_recovery.cc "PASSED")
endif()