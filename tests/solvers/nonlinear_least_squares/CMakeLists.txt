
##########
### LM ###
##########
if(PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_exe_and_test(lm_eigen solvers lm.cc "PASSED")
endif()
if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_serial_exe_and_test(lm_kokkos solvers lm_kokkos.cc "PASSED")
endif()

##########
### gn ###
##########
if(PRESSIO_ENABLE_TPL_EIGEN)
  set(fn gn_normal_eq_exp_data_fit_n2_eigen)
  add_serial_exe_and_test(${fn} solvers ${CMAKE_CURRENT_SOURCE_DIR}/${fn}.cc "PASSED")

  set(fn gn_normal_eq_hes_gra_api_rosenbrock4_eigen)
  add_serial_exe_and_test(${fn} solvers ${CMAKE_CURRENT_SOURCE_DIR}/${fn}.cc "PASSED")

  set(fn gn_normal_eq_res_jac_api_rosenbrock4_eigen)
  add_serial_exe_and_test(${fn} solvers ${CMAKE_CURRENT_SOURCE_DIR}/${fn}.cc "PASSED")

  set(fn gn_qr_line_search_rosenbrock3_eigen)
  add_serial_exe_and_test(${fn} solvers ${CMAKE_CURRENT_SOURCE_DIR}/${fn}.cc "PASSED")

  set(fn gn_qr_exp_data_fit_n2_eigen)
  add_serial_exe_and_test(${fn} solvers ${CMAKE_CURRENT_SOURCE_DIR}/${fn}.cc "PASSED")

  set(fn gn_qr_rosenbrock_scaled_eigen)
  add_serial_exe_and_test(${fn} solvers ${CMAKE_CURRENT_SOURCE_DIR}/${fn}.cc "PASSED")

  set(fn gn_qr_nist_mgh10)
  add_serial_exe_and_test(${fn} solvers ${CMAKE_CURRENT_SOURCE_DIR}/${fn}.cc "PASSED")
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
  set(f1 gn_normal_eq_exp_data_fit_n5_mpi_epetra)
  add_mpi_exe_and_test(${f1}_n2 solvers ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 2 "PASSED")
  add_mpi_exe_and_test(${f1}_n3 solvers ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 3 "PASSED")

  set(f1 gn_normal_eq_exp_data_fit_n11_mpi_epetra)
  add_mpi_exe_and_test(${f1}_n2 solvers ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 2 "PASSED")
  add_mpi_exe_and_test(${f1}_n3 solvers ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 3 "PASSED")

  set(f1 gn_qr_exp_data_fit_n5_mpi_epetra)
  add_mpi_exe_and_test(${f1}_n2 solvers ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 2 "PASSED")
  add_mpi_exe_and_test(${f1}_n3 solvers ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 3 "PASSED")

  set(f1 gn_qr_exp_data_fit_n11_mpi_epetra)
  add_mpi_exe_and_test(${f1}_n2 solvers ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 2 "PASSED")
  add_mpi_exe_and_test(${f1}_n3 solvers ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 3 "PASSED")
endif()
