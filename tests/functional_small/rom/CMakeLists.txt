
add_serial_utest(${TESTING_LEVEL}_rom_concepts ${CMAKE_CURRENT_SOURCE_DIR}/rom_concepts.cc)

if(PRESSIO_ENABLE_TPL_EIGEN)
 add_serial_utest(${TESTING_LEVEL}_rom_linear_decoder_arbitrary_types
 		${CMAKE_CURRENT_SOURCE_DIR}/linear_decoder_arbitrary_types.cc)

 add_serial_utest(${TESTING_LEVEL}_rom_linear_decoder_eigen
 		${CMAKE_CURRENT_SOURCE_DIR}/linear_decoder_eigen.cc)
endif()

add_subdirectory(correctness_galerkin)
add_subdirectory(correctness_lspg_steady)
add_subdirectory(correctness_lspg_unsteady)


if(PRESSIO_ENABLE_TPL_TRILINOS)
 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/hypred_updater_tpetra.cc)
 add_utest_mpi(${TESTING_LEVEL}_rom_hypred_updater_trilinos_tpetra "${SRC1}" gTestMain_tpetra 3)

 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/hypred_updater_epetra.cc)
 add_utest_mpi(${TESTING_LEVEL}_rom_hypred_updater_trilinos_epetra "${SRC1}" gTestMain_mpi 3)

 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/hypred_updater_tpetra_block.cc)
 add_utest_mpi(${TESTING_LEVEL}_rom_hypred_updater_trilinos_tpetra_block "${SRC1}" gTestMain_tpetra 3)
endif()
