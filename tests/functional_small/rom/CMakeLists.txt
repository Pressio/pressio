
add_serial_utest(${TESTING_LEVEL}_rom_concepts_A ${CMAKE_CURRENT_SOURCE_DIR}/concepts_A.cc)
add_serial_utest(${TESTING_LEVEL}_rom_concepts_B ${CMAKE_CURRENT_SOURCE_DIR}/concepts_B.cc)
#add_serial_utest(${TESTING_LEVEL}_rom_concepts_C ${CMAKE_CURRENT_SOURCE_DIR}/concepts_C.cc)
add_serial_utest(${TESTING_LEVEL}_rom_concepts_D ${CMAKE_CURRENT_SOURCE_DIR}/concepts_D.cc)

if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC linear_subspace_test1)
  add_serial_utest(${TESTING_LEVEL}_rom_${SRC} ${CMAKE_CURRENT_SOURCE_DIR}/${SRC}.cc)

  set(SRC trial_subspace_stdvec_full_state_eigen_reduced_state)
  add_serial_utest(${TESTING_LEVEL}_rom_${SRC} ${CMAKE_CURRENT_SOURCE_DIR}/${SRC}.cc)
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
 set(SRC trial_subspace_tpetra_full_state_eigen_reduced_state)
 add_utest_mpi(${TESTING_LEVEL}_rom_${SRC} "${SRC}.cc" gTestMain_tpetra 3)
endif()

add_subdirectory(galerkin_steady)
add_subdirectory(galerkin_unsteady)
add_subdirectory(lspg_steady)
add_subdirectory(lspg_unsteady)
