
if (PRESSIO_ENABLE_TPL_EIGEN)
add_serial_exe_and_test(newton_raphson_eigen solvers_nonlinear eigen.cc "PASSED")
endif()
