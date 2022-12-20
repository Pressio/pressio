
if(PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_utest(${TESTING_LEVEL}_rom_galerkin_unsteady_implicit_exe1
    ${CMAKE_CURRENT_SOURCE_DIR}/main1.cc)

  add_serial_utest(${TESTING_LEVEL}_rom_galerkin_unsteady_implicit_exe2
    ${CMAKE_CURRENT_SOURCE_DIR}/main2.cc)

  add_serial_utest(${TESTING_LEVEL}_rom_galerkin_unsteady_implicit_exe3
    ${CMAKE_CURRENT_SOURCE_DIR}/main3.cc)

  add_serial_utest(${TESTING_LEVEL}_rom_galerkin_unsteady_implicit_exe4
    ${CMAKE_CURRENT_SOURCE_DIR}/main4.cc)

  add_serial_utest(${TESTING_LEVEL}_rom_galerkin_unsteady_implicit_exe5
    ${CMAKE_CURRENT_SOURCE_DIR}/main5.cc)

  add_serial_utest(${TESTING_LEVEL}_rom_galerkin_unsteady_implicit_exe6
    ${CMAKE_CURRENT_SOURCE_DIR}/main6.cc)
endif()
