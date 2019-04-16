
IF(HAVE_TRILINOS)
  add_unit_test(sparse_matrix_epetra_class ${GTESTMAINSDIR}/gTestMain_mpi 3)
  add_unit_test(sparse_matrix_epetra_product_epetra_multi_vector ${GTESTMAINSDIR}/gTestMain_mpi 3)
  add_unit_test(sparse_matrix_tpetra_class ${GTESTMAINSDIR}/gTestMain_tpetra 3)
  add_unit_test(sparse_matrix_tpetra_traits ${GTESTMAINSDIR}/gTestMain_tpetra 3)
ENDIF()
