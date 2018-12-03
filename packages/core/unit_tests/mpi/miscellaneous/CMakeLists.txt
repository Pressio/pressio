
IF(HAVE_TRILINOS)
  add_unit_test(tpetra_maps_import ${GTESTMAINSDIR}/gTestMain_mpi 3)
ENDIF()