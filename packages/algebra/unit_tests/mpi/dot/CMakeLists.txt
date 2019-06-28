
IF(HAVE_TRILINOS)
  add_unit_test(multi_vector_dot_vector_distributed_epetra ${GTESTMAINSDIR}/gTestMain_mpi 3)
  add_unit_test(multi_vector_dot_multi_vector_distributed_epetra ${GTESTMAINSDIR}/gTestMain_mpi 3)

  add_unit_test(multi_vector_dot_vector_distributed_tpetra ${GTESTMAINSDIR}/gTestMain_tpetra 3)
  add_unit_test(multi_vector_dot_multi_vector_distributed_tpetra ${GTESTMAINSDIR}/gTestMain_tpetra 3)

ENDIF()
