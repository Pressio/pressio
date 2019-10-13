
# this function creates and adds an individual SERIAL test
macro(add_serial_exe_and_test TESTNAME PKGNAME TESTSRCS REGEX)
  # set name of the executable
  set(exeName ${PKGNAME}_${TESTNAME})
  add_executable(${exeName} ${TESTSRCS})
  add_test(NAME ${exeName} COMMAND ${exeName})
  set_tests_properties(
    ${exeName}
    PROPERTIES
    PASS_REGULAR_EXPRESSION ${REGEX}
  )

endmacro()
#=====================================================================

macro(add_mpi_exe_and_test TESTNAME PKGNAME TESTSRCS nRANKS REGEX)
  set(exeName ${PKGNAME}_${TESTNAME}_np${nRANKS})
  add_executable(${exeName} ${TESTSRCS})
  target_link_libraries(${exeName} ${MPI_CXX_LIBRARIES})
  add_test(
    NAME ${exeName}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nRANKS}
    ${MPIEXEC_PREFLAGS} ${exeName} ${MPIEXEC_POSTFLAGS}
    )
  set_tests_properties(
    ${exeName}
    PROPERTIES
    PASS_REGULAR_EXPRESSION ${REGEX}
  )
endmacro()
#=====================================================================
