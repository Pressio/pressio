
set(SRC1
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers_gn_test_conv_conditions.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers_gn_normal_eq_exp_data_fit_n2.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers_gn_qr_exp_data_fit_n2.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers_gn_qr_line_search_exp_data_fit_n2.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers_gn_qr_line_search_rosenbrock3.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers_gn_qr_line_search_rosenbrock4.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers_gn_normal_eq_rosenbrock4.cc)
add_serial_utest(solvers_nonlin_lsqr_serial "${SRC1}")
