
#add_serial_exe_and_test(exe1 ${TESTING_LEVEL}_rom_galerkin main1.cc "PASSED")

if(PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_utest(${TESTING_LEVEL}_rom_galerkin_exe1 ${CMAKE_CURRENT_SOURCE_DIR}/main1.cc)
  # add_serial_utest(${TESTING_LEVEL}_rom_galerkin_exe2 ${CMAKE_CURRENT_SOURCE_DIR}/main2.cc)
  # add_serial_utest(${TESTING_LEVEL}_rom_galerkin_exe3 ${CMAKE_CURRENT_SOURCE_DIR}/main3.cc)
  # add_serial_utest(${TESTING_LEVEL}_rom_galerkin_exe4 ${CMAKE_CURRENT_SOURCE_DIR}/main4.cc)

 # set(FN galerkin_default_with_mm_correctness_eigen)
 # add_serial_utest(${TESTING_LEVEL}_rom_${FN} ${CMAKE_CURRENT_SOURCE_DIR}/${FN}.cc)

 # # this is here because uses Eigen for rom types
 # set(FN galerkin_default_correctness_custom_types)
 # add_serial_utest(${TESTING_LEVEL}_rom_${FN} ${CMAKE_CURRENT_SOURCE_DIR}/${FN}.cc)

 # # this is here because uses Eigen for rom types
 # if(PRESSIO_ENABLE_TPL_TRILINOS)
 #  set(FN galerkin_default_correctness_tpetra)
 #  add_utest_mpi(${TESTING_LEVEL}_${FN} ${FN}.cc gTestMain_tpetra 3)
 # endif()
endif()
