
if(PRESSIO_ENABLE_TPL_TRILINOS)

  set(SRC1
    ${CMAKE_CURRENT_SOURCE_DIR}/tpetra/ops_multi_vector_distributed_tpetra_class.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/tpetra/ops_multi_vector_dot_multi_vector_distributed_tpetra.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/tpetra/ops_multi_vector_dot_self_distributed_tpetra.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/tpetra/ops_multi_vector_dot_vector_distributed_tpetra.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/tpetra/ops_multi_vector_tpetra_product_eigen_subview_vector_from_mv.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/tpetra/ops_multi_vector_tpetra_product_eigen_vector.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/tpetra/ops_multi_vector_tpetra_product_kokkos_vector.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/tpetra/ops_vector_distributed_tpetra_class.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/tpetra/ops_vector_distributed_tpetra_norm2.cc
    )
  add_utest_mpi(ops_tpetra "${SRC1}" gTestMain_tpetra 3 )

  set(SRC2
    ${CMAKE_CURRENT_SOURCE_DIR}/epetra/ops_multi_vector_epetra_product_eigen_vector.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/epetra/ops_vector_distributed_epetra_norm2.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/epetra/ops_multi_vector_dot_self_distributed_epetra.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/epetra/ops_multi_vector_dot_vector_distributed_epetra.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/epetra/ops_multi_vector_dot_multi_vector_distributed_epetra.cc
    )
  add_utest_mpi(ops_epetra "${SRC2}" gTestMain_mpi 3 )

endif()
