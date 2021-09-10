
if(PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_utest(${TESTING_LEVEL}_qr_eigen qr_eigen.cc)
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
  add_utest_mpi(${TESTING_LEVEL}_qr_epetra qr_epetra.cc gTestMain_mpi 2 )
  add_utest_mpi(${TESTING_LEVEL}_qr_tpetra qr_tpetra.cc gTestMain_tpetra 2 )
endif()