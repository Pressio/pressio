
IF(HAVE_TRILINOS)
  add_unit_test(multi_vector_distributed_epetra_factorization ${GTESTMAINSDIR}/gTestMain_mpi 2)
  add_unit_test(multi_vector_distributed_epetra_project ${GTESTMAINSDIR}/gTestMain_mpi 2)
  add_unit_test(multi_vector_distributed_epetra_solve ${GTESTMAINSDIR}/gTestMain_mpi 2)

  add_unit_test(multi_vector_distributed_tpetra_factorization ${GTESTMAINSDIR}/gTestMain_mpi 2)
  add_unit_test(multi_vector_distributed_tpetra_project ${GTESTMAINSDIR}/gTestMain_mpi 2)
  add_unit_test(multi_vector_distributed_tpetra_solve ${GTESTMAINSDIR}/gTestMain_mpi 2)
ENDIF()
