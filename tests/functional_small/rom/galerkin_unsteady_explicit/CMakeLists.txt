
if(PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_utest(${TESTING_LEVEL}_rom_galerkin_unsteady_explicit_exe1
    ${CMAKE_CURRENT_SOURCE_DIR}/main1.cc)
  add_serial_utest(${TESTING_LEVEL}_rom_galerkin_unsteady_explicit_exe2
    ${CMAKE_CURRENT_SOURCE_DIR}/main2.cc)
  add_serial_utest(${TESTING_LEVEL}_rom_galerkin_unsteady_explicit_exe3
    ${CMAKE_CURRENT_SOURCE_DIR}/main3.cc)
  add_serial_utest(${TESTING_LEVEL}_rom_galerkin_unsteady_explicit_exe5
    ${CMAKE_CURRENT_SOURCE_DIR}/main5.cc)
endif()
