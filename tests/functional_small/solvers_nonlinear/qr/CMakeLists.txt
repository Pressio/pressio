
if(PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_utest(${TESTING_LEVEL}_qr_eigen qr_eigen.cc)
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
  add_utest_mpi(${TESTING_LEVEL}_qr_tpetra gTestMain_tpetra 2 qr_tpetra.cc)
endif()
