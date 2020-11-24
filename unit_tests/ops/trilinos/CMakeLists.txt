
set(SRC1
  ${CMAKE_CURRENT_SOURCE_DIR}/tpetra/ops_vector_distributed_tpetra_norm2.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/tpetra/ops_vector_distributed_tpetra_setter.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/tpetra/ops_vector_distributed_tpetra_dot.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/tpetra/ops_asdiagonalmatrix_product_multivector.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/tpetra/ops_vector_elementwise_multiply.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/tpetra/ops_vector_pow.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/tpetra/ops_vector_abs.cc
)
add_utest_mpi(ops_tpetra "${SRC1}" gTestMain_tpetra 3)

set(SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/tpetra_block/ops_vector_pow.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/tpetra_block/ops_vector_elementwise_multiply.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/tpetra_block/ops_vector_abs.cc
)
add_utest_mpi(ops_tpetra_block "${SRC}" gTestMain_tpetra 3)

set(SRC2
  ${CMAKE_CURRENT_SOURCE_DIR}/epetra/ops_vector_distributed_epetra_dot.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/epetra/ops_vector_distributed_epetra_norm2.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/epetra/ops_multi_vector_epetra_product_eigen_vector.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/epetra/ops_vector_elementwise_multiply.cc
  )
add_utest_mpi(ops_epetra "${SRC2}" gTestMain_mpi 3 )
