
if(PRESSIO_ENABLE_INTERNAL_SPDLOG)

add_serial_utest(${TESTING_LEVEL}_logger logger.cc)

if(PRESSIO_ENABLE_TPL_MPI)
  add_utest_mpi(${TESTING_LEVEL}_logger_mpi gTestMain_mpi 2 logger_mpi.cc)
endif()

endif()

