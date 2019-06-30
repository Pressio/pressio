
IF(HAVE_TRILINOS)
  add_unit_test(vector_distributed_epetra_traits ${GTESTMAINSDIR}/gTestMain_mpi 1)
  add_unit_test(vector_distributed_epetra_class ${GTESTMAINSDIR}/gTestMain_mpi 3)
  add_unit_test(vector_distributed_tpetra_traits ${GTESTMAINSDIR}/gTestMain_tpetra 3)
  add_unit_test(vector_distributed_tpetra_class ${GTESTMAINSDIR}/gTestMain_tpetra 3)
  add_unit_test(vector_distributed_tpetra_block_traits ${GTESTMAINSDIR}/gTestMain_tpetra 3)
  add_unit_test(vector_distributed_tpetra_block_class ${GTESTMAINSDIR}/gTestMain_tpetra 3)
ENDIF()
