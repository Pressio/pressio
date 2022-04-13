
set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/odrom_tpetra.cc)
add_utest_mpi(${TESTING_LEVEL}_odrom_trilinos_tpetra "${SRC1}" gTestMain_tpetra 5)
