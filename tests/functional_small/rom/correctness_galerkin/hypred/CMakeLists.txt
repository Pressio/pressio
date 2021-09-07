
if(PRESSIO_ENABLE_TPL_EIGEN)
 set(FN galerkin_hypred_velo_correctness_eigen)
 add_serial_utest(${TESTING_LEVEL}_rom_${FN} ${CMAKE_CURRENT_SOURCE_DIR}/${FN}.cc)

 # this is here because uses Eigen for rom types
 set(FN galerkin_hypred_velo_correctness_custom_types)
 add_serial_utest(${TESTING_LEVEL}_rom_${FN} ${CMAKE_CURRENT_SOURCE_DIR}/${FN}.cc)
endif()
