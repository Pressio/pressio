
if(PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_exe_and_test(lm_eigen ${ROOTNAME} lm_eigen.cc "PASSED")
endif()
if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_serial_exe_and_test(lm_kokkos ${ROOTNAME} lm_kokkos.cc "PASSED")
endif()
