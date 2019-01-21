
IF(HAVE_TRILINOS)
  add_unit_test(gn_normal_eq_exp_data_fit_n5_mpi_epetra ${GTESTMAINSDIR}/gTestMain_mpi 2)
  add_unit_test(gn_qr_exp_data_fit_n5_mpi_epetra ${GTESTMAINSDIR}/gTestMain_mpi 2)
ENDIF()
