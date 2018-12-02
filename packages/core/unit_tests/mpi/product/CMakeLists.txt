
IF(HAVE_TRILINOS)
  add_unit_test(multi_vector_epetra_product_eigen_vector ${GTESTMAINSDIR}/gTestMain_mpi 3)
  add_unit_test(multi_vector_tpetra_product_eigen_vector ${GTESTMAINSDIR}/gTestMain_tpetra 3)

IF(HAVE_ARMADILLO)
  add_unit_test(multi_vector_epetra_product_armadillo_vector ${GTESTMAINSDIR}/gTestMain_mpi 3)
ENDIF()

ENDIF()
