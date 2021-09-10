
set(SRC
  ${CMAKE_CURRENT_SOURCE_DIR}/tpetra_block/ops_vector_pow.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/tpetra_block/ops_vector_elementwise_multiply.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/tpetra_block/ops_vector_abs.cc
)
add_utest_mpi(ops_tpetra_block "${SRC}" gTestMain_tpetra 3)
