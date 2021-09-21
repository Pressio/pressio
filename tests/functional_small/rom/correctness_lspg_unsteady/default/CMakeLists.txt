
if(PRESSIO_ENABLE_TPL_EIGEN)
 set(FN lspg_unsteady_default_correctness_eigen)
 add_serial_utest(${TESTING_LEVEL}_rom_${FN} ${CMAKE_CURRENT_SOURCE_DIR}/main_eigen.cc)

 # this is here because uses Eigen for rom types
 set(FN lspg_unsteady_default_correctness_custom_types)
 add_serial_utest(${TESTING_LEVEL}_rom_${FN} ${CMAKE_CURRENT_SOURCE_DIR}/main_custom_types.cc)
endif()
