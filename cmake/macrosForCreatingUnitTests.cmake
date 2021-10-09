
# this function creates and adds an individual SERIAL unit-test
macro(add_serial_utest TESTNAME)
  # set name of the executable
  set(testNameFinal ${TESTNAME})

  add_executable(${testNameFinal}
    ${ARGN} ${GTESTMAINSDIR}/gTestMain_serial.cc)

  find_library(GMOCKLIB NAMES gmock PATHS ${GTEST_ROOT}/lib/ NO_DEFAULT_PATH) 
  target_link_libraries(${testNameFinal}
    GTest::gtest ${GMOCKLIB} GTest::Main)

  add_test(NAME ${testNameFinal} COMMAND ${testNameFinal})
endmacro()
#=====================================================================

macro(add_serial_utest_kokkos TESTNAME TESTSRCS)
  # set name of the executable
  set(testNameFinal ${TESTNAME})

  add_executable(${testNameFinal}
    ${TESTSRCS} ${GTESTMAINSDIR}/gTestMain_kokkos.cc)

  target_link_libraries(${testNameFinal}
    GTest::GTest GTest::Main)

  add_test(NAME ${testNameFinal} COMMAND ${testNameFinal})
endmacro()
#=====================================================================


macro(add_utest_mpi TESTNAME TESTSRCS gMAIN nRANKS)
  set(testNameFinal ${TESTNAME}_np${nRANKS})

  add_executable(${testNameFinal}
    ${TESTSRCS} ${GTESTMAINSDIR}/${gMAIN}.cc)

  target_link_libraries(${testNameFinal}
    ${MPI_CXX_LIBRARIES}
    GTest::GTest GTest::Main)

  add_test(
    NAME ${testNameFinal}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nRANKS}
    ${MPIEXEC_PREFLAGS} ${testNameFinal} ${MPIEXEC_POSTFLAGS}
    )
endmacro()
#=====================================================================
