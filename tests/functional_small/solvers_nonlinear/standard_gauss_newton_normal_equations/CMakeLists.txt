
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../problems)

if(PRESSIO_ENABLE_TPL_EIGEN)
  # this test is here because uses Eigen vector/matrix for least-squares operators
  # while using custom types for the system residual/jacobian
  set(fname gn_normal_eq_res_jac_api_rosenbrock4_custom_types)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")

  set(fname gn_normal_eq_exp_data_fit_n2_eigen)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")

  set(fname gn_normal_eq_hes_gra_api_rosenbrock4_eigen)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")

  set(fname gn_normal_eq_res_jac_api_rosenbrock4_eigen)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")

  set(fname gn_normal_eq_custom_update_eigen)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")

  set(fname gn_normal_eq_observer_eigen)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")

  set(fname gn_normal_eq_nan_residual_eigen)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")

  set(fname gn_normal_eq_print_metrics_eigen)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
  set(f1 gn_normal_eq_exp_data_fit_n5_epetra)
  add_mpi_exe_and_test(${f1} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 2 "PASSED")
  add_mpi_exe_and_test(${f1} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 3 "PASSED")

  set(f1 gn_normal_eq_exp_data_fit_n11_epetra)
  add_mpi_exe_and_test(${f1} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 2 "PASSED")
  add_mpi_exe_and_test(${f1} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 3 "PASSED")
endif()