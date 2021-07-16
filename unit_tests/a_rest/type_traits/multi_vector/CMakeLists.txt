
set(exeRootName containers_multi_vector)

#############
### EIGEN ###
#############
if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/${exeRootName}_eigen_dynamic_class.cc)
  add_serial_utest(${exeRootName}_eigen_dynamic "${SRC1}")
endif()

##############
### KOKKOS ###
##############
if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_serial_utest_kokkos(
    ${exeRootName}_sharedmem_kokkos_mv
    ${CMAKE_CURRENT_SOURCE_DIR}/${exeRootName}_sharedmem_kokkos_mv.cc )
endif()

################
### TRILINOS ###
################
if(PRESSIO_ENABLE_TPL_TRILINOS)
  #### EPETRA ####
  add_serial_utest(
    ${exeRootName}_distributed_epetra_traits
    ${CMAKE_CURRENT_SOURCE_DIR}/${exeRootName}_distributed_epetra_traits.cc)

  set(SRC2 ${CMAKE_CURRENT_SOURCE_DIR}/${exeRootName}_distributed_epetra_class.cc)
  add_utest_mpi(${exeRootName}_epetra_class "${SRC2}" gTestMain_mpi 3 )

  #### TPETRA ####
  set(SRC1
    ${CMAKE_CURRENT_SOURCE_DIR}/${exeRootName}_distributed_tpetra_traits.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/${exeRootName}_distributed_tpetra_class.cc
    )
  add_utest_mpi(${exeRootName}_tpetra_class "${SRC1}" gTestMain_tpetra 3 )

  set(SRC1
    ${CMAKE_CURRENT_SOURCE_DIR}/${exeRootName}_distributed_tpetra_constructors_and_smf.cc)
  add_utest_mpi(${exeRootName}_tpetra_constr_and_smf "${SRC1}" gTestMain_tpetra 3 )

  #### TPETRA BLOCK ####
  set(SRC1
    ${CMAKE_CURRENT_SOURCE_DIR}/${exeRootName}_distributed_tpetra_block_traits.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/${exeRootName}_distributed_tpetra_block_class.cc)
  add_utest_mpi(${exeRootName}_tpetra_block_mv "${SRC1}" gTestMain_tpetra 3 )
endif()
