
IF(PRESSIO_ENABLE_TPL_TRILINOS)

  TRIBITS_ADD_EXECUTABLE(
    gn_normal_eq_exp_data_fit_n5_mpi_epetra
    SOURCES gn_normal_eq_exp_data_fit_n5_mpi_epetra.cc ${GTESTMAINSDIR}/gTestMain_mpi
    COMM mpi
    )
  TRIBITS_ADD_TEST(
    gn_normal_eq_exp_data_fit_n5_mpi_epetra
    COMM mpi  NUM_MPI_PROCS 2
    PASS_REGULAR_EXPRESSION "PASSED"  FAIL_REGULAR_EXPRESSION "FAILED"
    )
  TRIBITS_ADD_TEST(
    gn_normal_eq_exp_data_fit_n5_mpi_epetra
    COMM mpi  NUM_MPI_PROCS 3
    PASS_REGULAR_EXPRESSION "PASSED"  FAIL_REGULAR_EXPRESSION "FAILED"
    )
  #--------------------------------------------------------------


  TRIBITS_ADD_EXECUTABLE(
    gn_normal_eq_exp_data_fit_n11_mpi_epetra
    SOURCES gn_normal_eq_exp_data_fit_n11_mpi_epetra.cc ${GTESTMAINSDIR}/gTestMain_mpi
    COMM mpi
    )
  TRIBITS_ADD_TEST(
    gn_normal_eq_exp_data_fit_n11_mpi_epetra
    COMM mpi  NUM_MPI_PROCS 2
    PASS_REGULAR_EXPRESSION "PASSED"  FAIL_REGULAR_EXPRESSION "FAILED"
    )
  TRIBITS_ADD_TEST(
    gn_normal_eq_exp_data_fit_n11_mpi_epetra
    COMM mpi  NUM_MPI_PROCS 3
    PASS_REGULAR_EXPRESSION "PASSED"  FAIL_REGULAR_EXPRESSION "FAILED"
    )
  #--------------------------------------------------------------


  TRIBITS_ADD_EXECUTABLE(
    gn_qr_exp_data_fit_n5_mpi_epetra
    SOURCES gn_qr_exp_data_fit_n5_mpi_epetra.cc ${GTESTMAINSDIR}/gTestMain_mpi
    COMM mpi
    )
  TRIBITS_ADD_TEST(
    gn_qr_exp_data_fit_n5_mpi_epetra
    COMM mpi  NUM_MPI_PROCS 2
    PASS_REGULAR_EXPRESSION "PASSED"  FAIL_REGULAR_EXPRESSION "FAILED"
    )
  TRIBITS_ADD_TEST(
    gn_qr_exp_data_fit_n5_mpi_epetra
    COMM mpi  NUM_MPI_PROCS 3
    PASS_REGULAR_EXPRESSION "PASSED"  FAIL_REGULAR_EXPRESSION "FAILED"
    )
  #--------------------------------------------------------------


  TRIBITS_ADD_EXECUTABLE(
    gn_qr_exp_data_fit_n11_mpi_epetra
    SOURCES gn_qr_exp_data_fit_n11_mpi_epetra.cc ${GTESTMAINSDIR}/gTestMain_mpi
    COMM mpi
    )
  TRIBITS_ADD_TEST(
    gn_qr_exp_data_fit_n11_mpi_epetra
    COMM mpi  NUM_MPI_PROCS 2
    PASS_REGULAR_EXPRESSION "PASSED"  FAIL_REGULAR_EXPRESSION "FAILED"
    )
  TRIBITS_ADD_TEST(
    gn_qr_exp_data_fit_n11_mpi_epetra
    COMM mpi  NUM_MPI_PROCS 3
    PASS_REGULAR_EXPRESSION "PASSED"  FAIL_REGULAR_EXPRESSION "FAILED"
    )

ENDIF()
