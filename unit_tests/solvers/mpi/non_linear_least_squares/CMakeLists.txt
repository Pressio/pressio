
if(PRESSIO_ENABLE_TPL_TRILINOS)

  set(SRC1 solvers_gn_normal_eq_exp_data_fit_n5_mpi_epetra.cc)
  add_utest_mpi(solvers_gn_normal_eq_exp_data_fit_n5_mpi_epetra "${SRC1}" gTestMain_mpi 2 )
  add_utest_mpi(solvers_gn_normal_eq_exp_data_fit_n5_mpi_epetra "${SRC1}" gTestMain_mpi 3 )

  set(SRC2 solvers_gn_normal_eq_exp_data_fit_n11_mpi_epetra.cc)
  add_utest_mpi(solvers_gn_normal_eq_exp_data_fit_n11_mpi_epetra "${SRC2}" gTestMain_mpi 2 )
  add_utest_mpi(solvers_gn_normal_eq_exp_data_fit_n11_mpi_epetra "${SRC2}" gTestMain_mpi 3 )

  set(SRC3 solvers_gn_qr_exp_data_fit_n5_mpi_epetra.cc)
  add_utest_mpi(solvers_gn_qr_exp_data_fit_n5_mpi_epetra "${SRC3}" gTestMain_mpi 2 )
  add_utest_mpi(solvers_gn_qr_exp_data_fit_n5_mpi_epetra "${SRC3}" gTestMain_mpi 3 )

  set(SRC4 solvers_gn_qr_exp_data_fit_n11_mpi_epetra.cc)
  add_utest_mpi(solvers_gn_qr_exp_data_fit_n11_mpi_epetra "${SRC4}" gTestMain_mpi 2 )
  add_utest_mpi(solvers_gn_qr_exp_data_fit_n11_mpi_epetra "${SRC4}" gTestMain_mpi 3 )

endif()
