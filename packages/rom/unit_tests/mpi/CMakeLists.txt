
IF(HAVE_TRILINOS)
  add_unit_test(fom_states_epetra ${GTESTMAINSDIR}/gTestMain_mpi 2)
ENDIF()
