
If(HAVE_TRILINOS)
  add_unit_test(multi_vector_distributed_epetra_factorization_out_of_place ${GTESTMAINSDIR}/gTestMain_mpi 2)
  add_unit_test(multi_vector_distributed_tpetra_factorization_out_of_place ${GTESTMAINSDIR}/gTestMain_mpi 2)
  add_unit_test(multi_vector_distributed_epetra_tsqr_out_of_place_and_solve ${GTESTMAINSDIR}/gTestMain_mpi 2)
  add_unit_test(multi_vector_distributed_tpetra_tsqr_out_of_place_and_solve ${GTESTMAINSDIR}/gTestMain_mpi 2)
ENDIF()
