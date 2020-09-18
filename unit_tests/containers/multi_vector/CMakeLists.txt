
set(exeRootName containers_multi_vector)

set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/containers_multi_vector_serial_eigen_dynamic_class.cc)
add_serial_utest(${exeRootName}_serial "${SRC1}")


if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_serial_utest_kokkos(${exeRootName}_sharedmem_kokkos_class
    ${CMAKE_CURRENT_SOURCE_DIR}/${exeRootName}_sharedmem_kokkos_class.cc )
endif()


if(PRESSIO_ENABLE_TPL_TRILINOS)
  set(SRC1
    ${CMAKE_CURRENT_SOURCE_DIR}/${exeRootName}_distributed_tpetra_traits.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/${exeRootName}_distributed_tpetra_class.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/${exeRootName}_distributed_tpetra_block_traits.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/${exeRootName}_distributed_tpetra_block_class.cc
    )
  add_utest_mpi(containers_tpetra "${SRC1}" gTestMain_tpetra 3 )

  add_serial_utest(${exeRootName}_distributed_epetra_traits
    ${CMAKE_CURRENT_SOURCE_DIR}/${exeRootName}_distributed_epetra_traits.cc)

  set(SRC2 ${CMAKE_CURRENT_SOURCE_DIR}/${exeRootName}_distributed_epetra_class.cc)
  add_utest_mpi(containers_epetra "${SRC2}" gTestMain_mpi 3 )
endif()
