
macro(add_serial_utest TESTNAME)
  set(testNameFinal ${TESTNAME})
  add_executable(${testNameFinal} ${ARGN} ${GTESTMAINSDIR}/gTestMain_serial.cc)
  target_link_libraries(${testNameFinal} pressio gtest_main)
  add_test(NAME ${testNameFinal} COMMAND ${testNameFinal})
endmacro()
#=====================================================================

macro(add_serial_utest_kokkos TESTNAME)
  set(testNameFinal ${TESTNAME})
  add_executable(${testNameFinal} ${ARGN} ${GTESTMAINSDIR}/gTestMain_kokkos.cc)
  target_link_libraries(${testNameFinal} ${KOKKOS_LIBS} pressio gtest_main)
  add_test(NAME ${testNameFinal} COMMAND ${testNameFinal})
endmacro()
#=====================================================================

macro(add_utest_mpi TESTNAME TESTSRCS gMAIN nRANKS)
  set(testNameFinal ${TESTNAME}_np${nRANKS})
  add_executable(${testNameFinal} ${TESTSRCS} ${GTESTMAINSDIR}/${gMAIN}.cc)
  target_link_libraries(${testNameFinal} ${MPI_CXX_LIBRARIES} pressio gtest_main)
  add_test(
    NAME ${testNameFinal}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nRANKS}
    ${MPIEXEC_PREFLAGS} ${testNameFinal} ${MPIEXEC_POSTFLAGS}
    )
endmacro()
#=====================================================================
