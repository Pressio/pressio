
IF(HAVE_TRILINOS)
  add_unit_test(vector_distributed_epetra_norm2 ${GTESTMAINSDIR}/gTestMain_mpi 3)
  add_unit_test(vector_distributed_tpetra_norm2 ${GTESTMAINSDIR}/gTestMain_tpetra 3)
ENDIF()
