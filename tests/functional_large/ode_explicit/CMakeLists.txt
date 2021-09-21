
set(ROOTNAME ${TESTING_LEVEL}_ode_explicit)

if(PRESSIO_ENABLE_TPL_EIGEN)
add_serial_exe_and_test(simple_analytical_eigen ${ROOTNAME} simple_analytical_eigen.cc "PASSED")
endif()
