
if(PRESSIO_ENABLE_TPL_TRILINOS)

  set(SRC1
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_distributed_tpetra_traits.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_distributed_tpetra_class.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_distributed_tpetra_block_traits.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_distributed_tpetra_block_class.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_tpetra_maps_import.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_multi_vector_distributed_tpetra_traits.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_multi_vector_distributed_tpetra_class.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_multi_vector_distributed_tpetra_block_traits.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_multi_vector_distributed_tpetra_block_class.cc
    )
  add_utest_mpi(containers_tpetra "${SRC1}" gTestMain_tpetra 3 )

  add_utest_mpi(containers_vector_manipulation_tpetra
    containers_vector_manipulation_tpetra.cc gTestMain_tpetra 4)

  add_serial_utest(containers_vector_distributed_epetra_traits
    containers_vector_distributed_epetra_traits.cc)
  add_serial_utest(containers_multi_vector_distributed_epetra_traits
    containers_multi_vector_distributed_epetra_traits.cc)

  set(SRC2
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_distributed_epetra_class.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_multi_vector_distributed_epetra_class.cc
    )
  add_utest_mpi(containers_epetra "${SRC2}" gTestMain_mpi 3 )

endif()
