
#
# concepts testing
#
set(SRC concept_linear_subspace)
add_serial_utest(${TESTING_LEVEL}_${SRC} ${CMAKE_CURRENT_SOURCE_DIR}/${SRC}.cc)

# set(SRC concept_semidiscrete_fom)
# add_serial_utest(${TESTING_LEVEL}_${SRC} ${CMAKE_CURRENT_SOURCE_DIR}/${SRC}.cc)

# set(SRC concept_semidiscrete_fom_with_jacobian_action)
# add_serial_utest(${TESTING_LEVEL}_${SRC} ${CMAKE_CURRENT_SOURCE_DIR}/${SRC}.cc)

# set(SRC concept_semidiscrete_fom_with_mass_matrix_action)
# add_serial_utest(${TESTING_LEVEL}_${SRC} ${CMAKE_CURRENT_SOURCE_DIR}/${SRC}.cc)

#
# subspaces
#
if(PRESSIO_ENABLE_TPL_EIGEN)
set(SRC concept_possibly_affine_trial_subspace)
add_serial_utest(${TESTING_LEVEL}_${SRC} ${CMAKE_CURRENT_SOURCE_DIR}/${SRC}.cc)
endif()

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

#
# roms
#
# if(PRESSIO_ENABLE_TPL_EIGEN)
#   set(SRC reduced_default_operators_traits)
#   add_serial_utest(${TESTING_LEVEL}_rom_${SRC} ${CMAKE_CURRENT_SOURCE_DIR}/${SRC}.cc)
# endif()

add_subdirectory(galerkin_steady)
add_subdirectory(galerkin_unsteady_explicit)
add_subdirectory(galerkin_unsteady_implicit)
add_subdirectory(lspg_steady)
add_subdirectory(lspg_unsteady)
#add_subdirectory(lspg_unsteady_fully_discrete_with_step_recovery)
