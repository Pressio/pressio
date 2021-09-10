
add_serial_utest(${TESTING_LEVEL}_utils_serial_printer utils_serial_printer.cc)
add_serial_utest(${TESTING_LEVEL}_logger logger.cc)

if(PRESSIO_ENABLE_TPL_MPI)
  add_utest_mpi(${TESTING_LEVEL}_logger_mpi logger_mpi.cc gTestMain_mpi 2)
endif()
