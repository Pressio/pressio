
if(PRESSIO_ENABLE_TPL_EIGEN)
 set(FN lspg_steady_default_correctness)
 add_serial_utest(${TESTING_LEVEL}_rom_${FN} ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)
endif()
