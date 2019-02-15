
IF(HAVE_TRILINOS)
  add_unit_test(multi_vector_dot_self_distributed_epetra ${GTESTMAINSDIR}/gTestMain_mpi 3)
  add_unit_test(multi_vector_dot_self_distributed_tpetra ${GTESTMAINSDIR}/gTestMain_mpi 3)
ENDIF()
