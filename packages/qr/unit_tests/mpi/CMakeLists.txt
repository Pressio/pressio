
IF(HAVE_TRILINOS AND HAVE_ANASAZI_TSQR)
  #add_unit_test(multi_vector_distributed_epetra_factorization_in_place ${GTESTMAINSDIR}/gTestMain_mpi 2)
  add_unit_test(multi_vector_distributed_epetra_factorization_out_of_place ${GTESTMAINSDIR}/gTestMain_mpi 2)
  add_unit_test(multi_vector_distributed_epetra_out_of_place_and_solve ${GTESTMAINSDIR}/gTestMain_mpi 2)

  #add_unit_test(multi_vector_distributed_tpetra_factorization_in_place ${GTESTMAINSDIR}/gTestMain_mpi 2)
  add_unit_test(multi_vector_distributed_tpetra_factorization_out_of_place ${GTESTMAINSDIR}/gTestMain_mpi 2)
  add_unit_test(multi_vector_distributed_tpetra_out_of_place_and_solve ${GTESTMAINSDIR}/gTestMain_mpi 2)

  #add_unit_test(multi_vector_distributed_epetra_project ${GTESTMAINSDIR}/gTestMain_mpi 2)
  #add_unit_test(multi_vector_distributed_epetra_solve ${GTESTMAINSDIR}/gTestMain_mpi 2)
  #add_unit_test(multi_vector_distributed_tpetra_factorization ${GTESTMAINSDIR}/gTestMain_mpi 2)
  #add_unit_test(multi_vector_distributed_tpetra_project ${GTESTMAINSDIR}/gTestMain_mpi 2)
  #add_unit_test(multi_vector_distributed_tpetra_solve ${GTESTMAINSDIR}/gTestMain_mpi 2)
ENDIF()
