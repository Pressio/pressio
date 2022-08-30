
if(PRESSIO_ENABLE_TPL_EIGEN)
 # this is here because uses Eigen for rom types
 set(FN lspg_unsteady_prec_hyperreduced_correctness)
 add_serial_utest(${TESTING_LEVEL}_rom_${FN} ${CMAKE_CURRENT_SOURCE_DIR}/main.cc)

 # this is here because uses Eigen for rom types
 if(PRESSIO_ENABLE_TPL_TRILINOS)
   set(FN lspg_unsteady_prec_hyperreduced_correctness_tpetra)
   add_utest_mpi(${TESTING_LEVEL}_${FN} main_tpetra.cc gTestMain_tpetra 3)
 endif()
endif()
