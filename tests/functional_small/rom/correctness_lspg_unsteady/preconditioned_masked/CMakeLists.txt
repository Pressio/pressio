
if(PRESSIO_ENABLE_TPL_EIGEN)
 # this is here because uses Eigen for rom types
 set(FN lspg_unsteady_preconditioned_masked_correctness)
 add_serial_utest(${TESTING_LEVEL}_rom_${FN} ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)
endif()
