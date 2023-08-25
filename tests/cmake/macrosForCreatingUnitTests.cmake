
# this function creates and adds an individual SERIAL unit-test
macro(add_serial_utest TESTNAME)
  # set name of the executable
  set(testNameFinal ${TESTNAME})

  add_executable(${testNameFinal}
    ${ARGN} ${GTESTMAINSDIR}/gTestMain_serial.cc)

  target_link_libraries(${testNameFinal}
    GTest::gtest GTest::Main)

  add_test(NAME ${testNameFinal} COMMAND ${testNameFinal})
endmacro()
#=====================================================================
