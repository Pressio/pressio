

set(SRC1
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers_res_jac_api.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers_hess_grad_api.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers_valid_lin_solver.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers_valid_qr_solver.cc)
add_serial_utest(solvers_apis_check "${SRC1}")
