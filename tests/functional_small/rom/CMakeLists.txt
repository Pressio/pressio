
add_serial_utest(${TESTING_LEVEL}_rom_concepts ${CMAKE_CURRENT_SOURCE_DIR}/rom_concepts.cc)

if(PRESSIO_ENABLE_TPL_EIGEN)
 add_serial_utest(${TESTING_LEVEL}_rom_linear_decoder_arbitrary_types 
 		${CMAKE_CURRENT_SOURCE_DIR}/linear_decoder_arbitrary_types.cc)

 add_serial_utest(${TESTING_LEVEL}_rom_linear_decoder_eigen 
 		${CMAKE_CURRENT_SOURCE_DIR}/linear_decoder_eigen.cc)
endif()

add_subdirectory(galerkin_default_correctness)
add_subdirectory(galerkin_masked_velocity_correctness)
add_subdirectory(galerkin_hypred_velocity_correctness)