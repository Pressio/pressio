
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../problems)

if(PRESSIO_ENABLE_TPL_EIGEN)
  set(fname gn_qr_line_search_rosenbrock3_eigen)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")

  set(fname gn_qr_exp_data_fit_n2_eigen)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")

  set(fname gn_qr_rosenbrock_scaled_eigen)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")

  set(fname gn_qr_nist_mgh10_eigen)
  add_serial_exe_and_test(${fname} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${fname}.cc "PASSED")
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
  set(f1 gn_qr_exp_data_fit_n5_epetra)
  add_mpi_exe_and_test(${f1} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 2 "PASSED")
  add_mpi_exe_and_test(${f1} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 3 "PASSED")

  set(f1 gn_qr_exp_data_fit_n11_epetra)
  add_mpi_exe_and_test(${f1} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 2 "PASSED")
  add_mpi_exe_and_test(${f1} ${ROOTNAME} ${CMAKE_CURRENT_SOURCE_DIR}/${f1}.cc 3 "PASSED")
endif()