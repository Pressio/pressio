
if (PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_exe_and_test(gn_irwls_correctness solvers_nonlinear main.cc "PASSED")
endif()
