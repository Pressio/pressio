
if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_serial_exe_and_test(geqrf_kokkos solvers geqrf_kokkos.cc "PASSED")
  add_serial_exe_and_test(getrs_kokkos solvers getrs_kokkos.cc "PASSED")
endif()

if(PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_exe_and_test(LSCG_iterative_dense_eigen solvers solvers_linear_dense_eigen.cc "PASSED")
  target_compile_definitions(solvers_LSCG_iterative_dense_eigen PRIVATE DO_LSCG)

  add_serial_exe_and_test(ColPivHSQR_direct_dense_eigen solvers solvers_linear_dense_eigen.cc "PASSED")
  target_compile_definitions(solvers_ColPivHSQR_direct_dense_eigen PRIVATE DO_ColPivHSQR)

  add_serial_exe_and_test(HouseholderQR_direct_dense_eigen solvers solvers_linear_dense_eigen.cc "PASSED")
  target_compile_definitions(solvers_HouseholderQR_direct_dense_eigen PRIVATE DO_HouseholderQR)

  add_serial_exe_and_test(LSCG_iterative_sparse_eigen solvers solvers_linear_sparse_eigen.cc "PASSED")
endif()
