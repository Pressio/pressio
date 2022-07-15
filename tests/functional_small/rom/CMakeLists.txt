#
# concepts
#
add_serial_utest(${TESTING_LEVEL}_rom_concepts_A ${CMAKE_CURRENT_SOURCE_DIR}/concepts_A.cc)
add_serial_utest(${TESTING_LEVEL}_rom_concepts_B ${CMAKE_CURRENT_SOURCE_DIR}/concepts_B.cc)
add_serial_utest(${TESTING_LEVEL}_rom_concepts_C ${CMAKE_CURRENT_SOURCE_DIR}/concepts_C.cc)

#
# linear subspace
#
if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC trial_subspace_stdvec_full_state_eigen_reduced_state)
  add_serial_utest(${TESTING_LEVEL}_rom_${SRC} ${CMAKE_CURRENT_SOURCE_DIR}/${SRC}.cc)
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
  set(SRC trial_subspace_tpetra_full_state_eigen_reduced_state)
  add_utest_mpi(${TESTING_LEVEL}_rom_${SRC} "${SRC}.cc" gTestMain_tpetra 3)
endif()

add_subdirectory(galerkin_steady)
add_subdirectory(lspg_steady)








# if(PRESSIO_ENABLE_TPL_EIGEN)
#  add_serial_utest(${TESTING_LEVEL}_rom_linear_decoder_arbitrary_types
#    ${CMAKE_CURRENT_SOURCE_DIR}/linear_decoder_arbitrary_types.cc)
#  add_serial_utest(${TESTING_LEVEL}_rom_linear_decoder_eigen
#    ${CMAKE_CURRENT_SOURCE_DIR}/linear_decoder_eigen.cc)
#  add_serial_utest(${TESTING_LEVEL}_rom_residual_evaluator_eigen_bdf1
#    ${CMAKE_CURRENT_SOURCE_DIR}/residual_evaluator_eigen_bdf1.cc)
#  add_serial_utest(${TESTING_LEVEL}_rom_residual_evaluator_eigen_bdf2
#    ${CMAKE_CURRENT_SOURCE_DIR}/residual_evaluator_eigen_bdf2.cc)
#  add_serial_utest(${TESTING_LEVEL}_rom_residual_evaluator_eigen_cn
#    ${CMAKE_CURRENT_SOURCE_DIR}/residual_evaluator_eigen_cn.cc)
# endif()

# if(PRESSIO_ENABLE_TPL_TRILINOS)
# add_subdirectory(correctness_lspg_steady)
# add_subdirectory(correctness_lspg_unsteady)
#  add_subdirectory(odrom_tpetra_constant_modes_per_tile)
#  add_subdirectory(odrom_tpetra_varying_modes_per_tile)
#  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/hypred_updater_tpetra.cc)
#  add_utest_mpi(${TESTING_LEVEL}_rom_hypred_updater_trilinos_tpetra "${SRC1}" gTestMain_tpetra 3)
#  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/hypred_updater_epetra.cc)
#  add_utest_mpi(${TESTING_LEVEL}_rom_hypred_updater_trilinos_epetra "${SRC1}" gTestMain_mpi 3)
#  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/hypred_updater_tpetra_block.cc)
#  add_utest_mpi(${TESTING_LEVEL}_rom_hypred_updater_trilinos_tpetra_block "${SRC1}" gTestMain_tpetra 3)
# endif()
