
IF(HAVE_TRILINOS)
  add_unit_test(vector_distributed_epetra_traits ${GTESTMAINSDIR}/gTestMain_mpi 1)
  add_unit_test(vector_distributed_epetra_class ${GTESTMAINSDIR}/gTestMain_mpi 3)
  add_unit_test(multi_vector_distributed_epetra_traits ${GTESTMAINSDIR}/gTestMain_mpi 1)
  add_unit_test(multi_vector_distributed_epetra_class ${GTESTMAINSDIR}/gTestMain_mpi 3)
  add_unit_test(multi_vector_dot_vector_distributed_epetra_class ${GTESTMAINSDIR}/gTestMain_mpi 3)
  add_unit_test(multi_vector_epetra_product_eigen_vector ${GTESTMAINSDIR}/gTestMain_mpi 3)
  add_unit_test(sparse_matrix_epetra_class ${GTESTMAINSDIR}/gTestMain_mpi 3)
  add_unit_test(sparse_matrix_epetra_product_epetra_multi_vector ${GTESTMAINSDIR}/gTestMain_mpi 3)

  add_unit_test(vector_distributed_tpetra_traits ${GTESTMAINSDIR}/gTestMain_tpetra 3)
  add_unit_test(vector_distributed_tpetra_class ${GTESTMAINSDIR}/gTestMain_tpetra 3)
  add_unit_test(multi_vector_distributed_tpetra_traits ${GTESTMAINSDIR}/gTestMain_tpetra 3)
  add_unit_test(multi_vector_distributed_tpetra_class ${GTESTMAINSDIR}/gTestMain_tpetra 3)
  add_unit_test(multi_vector_tpetra_product_eigen_vector ${GTESTMAINSDIR}/gTestMain_tpetra 3)
  add_unit_test(multi_vector_dot_vector_distributed_tpetra_class ${GTESTMAINSDIR}/gTestMain_mpi 3)
ENDIF()

IF(HAVE_TRILINOS AND HAVE_ARMADILLO)
  add_unit_test(multi_vector_epetra_product_armadillo_vector ${GTESTMAINSDIR}/gTestMain_mpi 3)
ENDIF()

