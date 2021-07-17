
if (PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_exe_and_test   (irwls_gn_correctness solvers main.cc "PASSED")
endif()
