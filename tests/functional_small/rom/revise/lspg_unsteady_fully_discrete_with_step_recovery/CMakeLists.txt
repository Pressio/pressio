
if(PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_utest(${TESTING_LEVEL}_rom_lspg_unsteady_fully_discrete_api_step_recovery_exe2states
    ${CMAKE_CURRENT_SOURCE_DIR}/main_2states.cc)

  add_serial_utest(${TESTING_LEVEL}_rom_lspg_unsteady_fully_discrete_api_step_recovery_exe3states
    ${CMAKE_CURRENT_SOURCE_DIR}/main_3states.cc)
endif()
