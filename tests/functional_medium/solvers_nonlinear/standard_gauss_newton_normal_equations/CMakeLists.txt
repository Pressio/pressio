
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../problems)

if(PRESSIO_ENABLE_TPL_EIGEN)
  set(fname eigen_gn_normal_eq_exp_data_fit_n2)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")

  set(fname eigen_gn_normal_eq_hes_gra_api_rosenbrock4)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")

  set(fname eigen_gn_normal_eq_res_jac_api_rosenbrock4)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")

  set(fname eigen_gn_normal_eq_custom_update)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")

  set(fname eigen_gn_normal_eq_observer)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")

  set(fname eigen_gn_normal_eq_nan_residual)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")

  set(fname eigen_gn_normal_eq_print_metrics)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
  set(f1 epetra_gn_normal_eq_exp_data_fit_n5)
  add_mpi_exe_and_test(${f1} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 2 "PASSED")
  add_mpi_exe_and_test(${f1} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 3 "PASSED")

  set(f1 epetra_gn_normal_eq_exp_data_fit_n11)
  add_mpi_exe_and_test(${f1} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 2 "PASSED")
  add_mpi_exe_and_test(${f1} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 3 "PASSED")
endif()