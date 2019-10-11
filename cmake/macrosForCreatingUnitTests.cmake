
# this function creates and adds an individual SERIAL unit-test
macro(add_serial_utest TESTNAME TESTSRCS)
  # set name of the executable
  set(exeName utest_${TESTNAME})

  add_executable(${exeName}
    ${TESTSRCS} ${GTESTMAINSDIR}/gTestMain_serial.cc)

  target_link_libraries(${exeName}
    GTest::GTest GTest::Main)

  add_test(NAME ${exeName} COMMAND ${exeName})
endmacro()
#=====================================================================


macro(add_utest_mpi TESTNAME TESTSRCS gMAIN nRANKS)
  set(exeName utest_${TESTNAME}_np${nRANKS})

  add_executable(${exeName}
    ${TESTSRCS} ${GTESTMAINSDIR}/${gMAIN}.cc)

  target_link_libraries(${exeName}
    ${MPI_CXX_LIBRARIES}
    GTest::GTest GTest::Main)

  add_test(
    NAME ${exeName}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nRANKS}
    ${MPIEXEC_PREFLAGS} ${exeName} ${MPIEXEC_POSTFLAGS}
    )
endmacro()
#=====================================================================


# function(add_multi_utest_mpi TESTNAME gMAIN nRANKS SRCS)
#   set(exeName utest_${TESTNAME}_np${nRANKS})

#   add_executable(${exeName}
#     ${SRCS} ${GTESTMAINSDIR}/${gMAIN}.cc)

#   target_link_libraries(${exeName}
#     ${MPI_CXX_LIBRARIES} GTest::GTest GTest::Main)

#   add_test(
#     NAME ${exeName}
#     COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${nRANKS}
#     ${MPIEXEC_PREFLAGS} ${exeName} ${MPIEXEC_POSTFLAGS}
#     )
# endfunction()
