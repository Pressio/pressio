
macro(add_serial_exe_and_test TESTNAME PKGNAME TESTSRCS REGEX)
  set(testNameFinal ${PKGNAME}_${TESTNAME})
  add_executable(${testNameFinal} ${TESTSRCS})
  target_include_directories(${testNameFinal} PRIVATE ${CMAKE_SOURCE_DIR}/include)
  add_test(NAME ${testNameFinal} COMMAND ${testNameFinal})
  set_tests_properties(
    ${testNameFinal}
    PROPERTIES
    PASS_REGULAR_EXPRESSION ${REGEX}
    FAIL_REGULAR_EXPRESSION "FAILED"
  )
endmacro()
#=====================================================================

macro(add_mpi_exe_and_test TESTNAME PKGNAME TESTSRCS nRANKS REGEXOK)
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
    PASS_REGULAR_EXPRESSION ${REGEXOK}
    FAIL_REGULAR_EXPRESSION "FAILED"
  )
endmacro()
#=====================================================================
