
if(PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_utest(${TESTING_LEVEL}_rom_lspg_unsteady_exe1
    ${CMAKE_CURRENT_SOURCE_DIR}/main1.cc)

  add_serial_utest(${TESTING_LEVEL}_rom_lspg_unsteady_exe2
    ${CMAKE_CURRENT_SOURCE_DIR}/main2.cc)

  add_serial_utest(${TESTING_LEVEL}_rom_lspg_unsteady_exe3
    ${CMAKE_CURRENT_SOURCE_DIR}/main3.cc)
endif()
