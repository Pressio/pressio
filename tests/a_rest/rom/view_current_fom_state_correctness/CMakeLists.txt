
add_serial_exe_and_test(galerkin_view_current_fom_state_continuous_time_api
  rom ${CMAKE_CURRENT_SOURCE_DIR}/galerkin_contiuous_time_api.cc "PASSED")

add_serial_exe_and_test(galerkin_view_current_fom_state_discrete_time_api
  rom ${CMAKE_CURRENT_SOURCE_DIR}/galerkin_discrete_time_api.cc "PASSED")

add_serial_exe_and_test(lspg_view_current_fom_state_steady_api
  rom ${CMAKE_CURRENT_SOURCE_DIR}/lspg_steady_api.cc "PASSED")

add_serial_exe_and_test(lspg_view_current_continuous_time_api
  rom ${CMAKE_CURRENT_SOURCE_DIR}/lspg_continuous_time_api.cc "PASSED")

add_serial_exe_and_test(lspg_view_current_discrete_time_api
  rom ${CMAKE_CURRENT_SOURCE_DIR}/lspg_discrete_time_api.cc "PASSED")
