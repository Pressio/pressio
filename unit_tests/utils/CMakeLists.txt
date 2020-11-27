
add_serial_utest(utils_serial_printer utils_serial_printer.cc)
add_serial_utest(logger logger.cc)

if(PRESSIO_ENABLE_TPL_MPI)
  add_utest_mpi(logger_mpi logger_mpi.cc gTestMain_mpi 2)
endif()
