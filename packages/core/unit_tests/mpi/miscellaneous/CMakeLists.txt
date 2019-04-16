
IF(HAVE_TRILINOS)
  add_unit_test(tpetra_maps_import ${GTESTMAINSDIR}/gTestMain_tpetra 3)
  add_unit_test(mpi_print ${GTESTMAINSDIR}/gTestMain_mpi 2)
  add_unit_test(vector_manipulation_tpetra ${GTESTMAINSDIR}/gTestMain_tpetra 4)
ENDIF()
