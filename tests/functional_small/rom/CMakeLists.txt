
add_serial_utest(${TESTING_LEVEL}_rom_concepts ${CMAKE_CURRENT_SOURCE_DIR}/rom_concepts.cc)

if(PRESSIO_ENABLE_TPL_EIGEN)
 add_serial_utest(${TESTING_LEVEL}_rom_linear_decoder_arbitrary_types 
 		${CMAKE_CURRENT_SOURCE_DIR}/linear_decoder_arbitrary_types.cc)

 add_serial_utest(${TESTING_LEVEL}_rom_linear_decoder_eigen 
 		${CMAKE_CURRENT_SOURCE_DIR}/linear_decoder_eigen.cc)
endif()

# if(PRESSIO_ENABLE_TPL_TRILINOS)
#   add_utest_mpi(${TESTING_LEVEL}_qr_epetra qr_epetra.cc gTestMain_mpi 2 )
#   add_utest_mpi(${TESTING_LEVEL}_qr_tpetra qr_tpetra.cc gTestMain_tpetra 2 )
# endif()