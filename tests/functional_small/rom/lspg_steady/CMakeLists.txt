
if(PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_utest(${TESTING_LEVEL}_rom_lspg_steady_exe1 ${CMAKE_CURRENT_SOURCE_DIR}/main1.cc)
  add_serial_utest(${TESTING_LEVEL}_rom_lspg_steady_exe2 ${CMAKE_CURRENT_SOURCE_DIR}/main2.cc)
endif()
