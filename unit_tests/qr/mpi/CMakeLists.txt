
if(PRESSIO_ENABLE_TPL_TRILINOS)

  set(SRC1
    ${CMAKE_CURRENT_SOURCE_DIR}/qr_multi_vector_distributed_epetra_factorization_out_of_place.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/qr_multi_vector_distributed_epetra_tsqr_out_of_place_and_solve.cc
    )
  add_utest_mpi(qr_epetra "${SRC1}" gTestMain_mpi 2 )

  set(SRC2
    ${CMAKE_CURRENT_SOURCE_DIR}/qr_multi_vector_distributed_tpetra_factorization_out_of_place.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/qr_multi_vector_distributed_tpetra_tsqr_out_of_place_and_solve.cc
    )
  add_utest_mpi(qr_tpetra "${SRC2}" gTestMain_tpetra 2 )

endif()
