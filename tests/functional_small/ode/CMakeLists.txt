
add_serial_utest(
	${TESTING_LEVEL}_ode_explicit_concepts 
	${CMAKE_CURRENT_SOURCE_DIR}/ode_explicit_concepts.cc)

# if(PRESSIO_ENABLE_TPL_EIGEN)
# set(SRC2 ${CMAKE_CURRENT_SOURCE_DIR}/ode_explicit_eigen.cc)
# add_serial_utest(${TESTING_LEVEL}_ode_explicit ${SRC2})
# endif()
