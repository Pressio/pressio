
# this function creates and adds an individual SERIAL unit-test
#
# DIRNAME is where the TESTNAME.cc file lives
#
macro(add_serial_utest TESTNAME DIRNAME)
	set(exeName utest_${TESTNAME})
	add_executable(${exeName} ${DIRNAME}/${TESTNAME}.cc ${GTESTMAINSDIR}/gTestMain_serial.cc)
    target_link_libraries(${exeName} GTest::GTest GTest::Main)
    add_test(NAME ${exeName} COMMAND ${exeName})
endmacro()
