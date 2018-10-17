

IF(HAVE_TRILINOS)
  add_unit_test(multi_vector_distributed_epetra_qr_factorization ${GTESTMAINSDIR}/gTestMain_mpi 3)
ENDIF()
