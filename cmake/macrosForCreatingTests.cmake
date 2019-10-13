
# this function creates and adds an individual SERIAL test
macro(add_serial_exe_and_test TESTNAME PKGNAME TESTSRCS REGEX)
  # set name of the executable
  set(testNameFinal ${PKGNAME}_${TESTNAME})
  add_executable(${testNameFinal} ${TESTSRCS})
  add_test(NAME ${testNameFinal} COMMAND ${testNameFinal})
  set_tests_properties(
    ${testNameFinal}
    PROPERTIES
    PASS_REGULAR_EXPRESSION ${REGEX}
  )

endmacro()
#=====================================================================

macro(add_mpi_exe_and_test TESTNAME PKGNAME TESTSRCS nRANKS REGEX)
  set(testNameFinal ${PKGNAME}_${TESTNAME}_np${nRANKS})
  add_executable(${testNameFinal} ${TESTSRCS})
  target_link_libraries(${testNameFinal} ${MPI_CXX_LIBRARIES})
  add_test(
    NAME ${testNameFinal}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nRANKS}
    ${MPIEXEC_PREFLAGS} ${testNameFinal} ${MPIEXEC_POSTFLAGS}
    )
  set_tests_properties(
    ${testNameFinal}
    PROPERTIES
    PASS_REGULAR_EXPRESSION ${REGEX}
  )
endmacro()
#=====================================================================
