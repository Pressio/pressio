
if(PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_exe_and_test(lm_eigen solvers_nonlinear lm_eigen.cc "PASSED")
endif()
if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_serial_exe_and_test(lm_kokkos solvers_nonlinear lm_kokkos.cc "PASSED")
endif()
