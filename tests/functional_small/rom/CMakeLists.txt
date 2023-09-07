
set(ROM_SOURCES_CONCEPTS ${CMAKE_CURRENT_SOURCE_DIR}/concept_linear_subspace.cc)
add_serial_utest(${TESTING_LEVEL}_rom_concepts ${ROM_SOURCES_CONCEPTS})

if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC concept_trial_subspace)
  add_serial_utest(${TESTING_LEVEL}_${SRC} ${CMAKE_CURRENT_SOURCE_DIR}/${SRC}.cc)
endif()

#
# subspaces
#
if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC linear_subspace_test1)
  add_serial_utest(${TESTING_LEVEL}_rom_${SRC} ${CMAKE_CURRENT_SOURCE_DIR}/${SRC}.cc)
  set(SRC trial_subspace_stdvec_full_state_eigen_reduced_state)
  add_serial_utest(${TESTING_LEVEL}_rom_${SRC} ${CMAKE_CURRENT_SOURCE_DIR}/${SRC}.cc)
endif()


if(PRESSIO_ENABLE_TPL_TRILINOS AND PRESSIO_ENABLE_TPL_EIGEN)
 set(SRC trial_subspace_tpetra_full_state_eigen_reduced_state)
 add_utest_mpi(${TESTING_LEVEL}_rom_${SRC} "${SRC}.cc" gTestMain_tpetra 3)
endif()

#
# roms
#
if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SOURCES_GALERKIN_STEADY
    ${CMAKE_CURRENT_SOURCE_DIR}/galerkin_steady/main1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/galerkin_steady/main3.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/galerkin_steady/main4.cc)
  add_serial_utest(${TESTING_LEVEL}_rom_galerkin_steady ${SOURCES_GALERKIN_STEADY})

  set(SOURCES_GALERKIN_UNSTEADY_EXP
    ${CMAKE_CURRENT_SOURCE_DIR}/galerkin_unsteady_explicit/main1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/galerkin_unsteady_explicit/main2.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/galerkin_unsteady_explicit/main3.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/galerkin_unsteady_explicit/main5.cc)
  add_serial_utest(${TESTING_LEVEL}_rom_galerkin_unsteady_explicit ${SOURCES_GALERKIN_UNSTEADY_EXP})

  set(SOURCES_GALERKIN_UNSTEADY_IMP
    ${CMAKE_CURRENT_SOURCE_DIR}/galerkin_unsteady_implicit/main1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/galerkin_unsteady_implicit/main2.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/galerkin_unsteady_implicit/main3.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/galerkin_unsteady_implicit/main4.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/galerkin_unsteady_implicit/main5.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/galerkin_unsteady_implicit/main6.cc)
  add_serial_utest(${TESTING_LEVEL}_rom_galerkin_unsteady_implicit ${SOURCES_GALERKIN_UNSTEADY_IMP})

  set(SOURCES_LSPG_STEADY
    ${CMAKE_CURRENT_SOURCE_DIR}/lspg_steady/main1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/lspg_steady/main2.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/lspg_steady/main3.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/lspg_steady/main4.cc)
  add_serial_utest(${TESTING_LEVEL}_rom_lspg_steady ${SOURCES_LSPG_STEADY})

  set(SOURCES_LSPG_UNSTEADY
    ${CMAKE_CURRENT_SOURCE_DIR}/lspg_unsteady/main1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/lspg_unsteady/main2.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/lspg_unsteady/main3.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/lspg_unsteady/main4.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/lspg_unsteady/main5.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/lspg_unsteady/main6.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/lspg_unsteady/main7.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/lspg_unsteady/main8.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/lspg_unsteady/main9.cc)
  add_serial_utest(${TESTING_LEVEL}_rom_lspg_unsteady ${SOURCES_LSPG_UNSTEADY})

  set(SOURCES_MIXED_FOM_LSPG_UNSTEADY
    ${CMAKE_CURRENT_SOURCE_DIR}/lspg_unsteady_mixed_fom/main1_bdf1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/lspg_unsteady_mixed_fom/main1_bdf2.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/lspg_unsteady_mixed_fom/main2_bdf1.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/lspg_unsteady_mixed_fom/main2_bdf2.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/lspg_unsteady_mixed_fom/main3_bdf2.cc)
  add_serial_utest(${TESTING_LEVEL}_rom_lspg_mixed_fom_unsteady ${SOURCES_MIXED_FOM_LSPG_UNSTEADY})

endif()
