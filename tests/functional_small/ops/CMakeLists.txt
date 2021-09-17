
# to include google test fixtures
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/fixtures)

if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_eigen_vector.cc)
  add_serial_utest(${TESTING_LEVEL}_ops_vector_eigen "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_eigen_matrix.cc)
  add_serial_utest(${TESTING_LEVEL}_ops_matrix_eigen "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_eigen_diag.cc)
  add_serial_utest(${TESTING_LEVEL}_ops_diag_eigen "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_eigen_span.cc)
  add_serial_utest(${TESTING_LEVEL}_ops_span_eigen "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_eigen_subspan.cc)
  add_serial_utest(${TESTING_LEVEL}_ops_subspan_eigen "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_eigen_level2.cc)
  add_serial_utest(${TESTING_LEVEL}_ops_level2_eigen "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_eigen_level3.cc)
  add_serial_utest(${TESTING_LEVEL}_ops_level3_eigen "${SRC1}")
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_kokkos_vector.cc)
  add_serial_utest_kokkos(${TESTING_LEVEL}_ops_vector_kokkos "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_kokkos_matrix.cc)
  add_serial_utest_kokkos(${TESTING_LEVEL}_ops_matrix_kokkos "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_kokkos_diag.cc)
  add_serial_utest_kokkos(${TESTING_LEVEL}_ops_diag_kokkos "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_kokkos_span.cc)
  add_serial_utest_kokkos(${TESTING_LEVEL}_ops_span_kokkos "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_kokkos_subspan.cc)
  add_serial_utest_kokkos(${TESTING_LEVEL}_ops_subspan_kokkos "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_kokkos_level2.cc)
  add_serial_utest_kokkos(${TESTING_LEVEL}_ops_level2_kokkos "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_kokkos_level3.cc)
  add_serial_utest_kokkos(${TESTING_LEVEL}_ops_level3_kokkos "${SRC1}")
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
 # TPETRA
 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_tpetra_vector.cc)
 add_utest_mpi(${TESTING_LEVEL}_ops_vector_tpetra "${SRC1}" gTestMain_tpetra 3)

 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_tpetra_multi_vector.cc)
 add_utest_mpi(${TESTING_LEVEL}_ops_multi_vector_tpetra "${SRC1}" gTestMain_tpetra 3)

 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_tpetra_level2.cc)
 add_utest_mpi(${TESTING_LEVEL}_ops_level2_tpetra "${SRC1}" gTestMain_tpetra 3)

 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_tpetra_level3.cc)
 add_utest_mpi(${TESTING_LEVEL}_ops_level3_tpetra "${SRC1}" gTestMain_tpetra 3)

 # TPETRA BLOCK
 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_tpetra_block_vector.cc)
 add_utest_mpi(${TESTING_LEVEL}_ops_vector_tpetra_block "${SRC1}" gTestMain_tpetra 3)

 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_tpetra_block_multi_vector.cc)
 add_utest_mpi(${TESTING_LEVEL}_ops_multi_vector_tpetra_block "${SRC1}" gTestMain_tpetra 3)

 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_tpetra_block_level2.cc)
 add_utest_mpi(${TESTING_LEVEL}_ops_level2_tpetra_block "${SRC1}" gTestMain_tpetra 3)

 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_tpetra_block_level3.cc)
 add_utest_mpi(${TESTING_LEVEL}_ops_level3_tpetra_block "${SRC1}" gTestMain_tpetra 3)

 # EPETRA
 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_epetra_vector.cc)
 add_utest_mpi(${TESTING_LEVEL}_ops_vector_epetra "${SRC1}" gTestMain_mpi 3)

 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_epetra_multi_vector.cc)
 add_utest_mpi(${TESTING_LEVEL}_ops_multi_vector_epetra "${SRC1}" gTestMain_mpi 3)

 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_epetra_level2.cc)
 add_utest_mpi(${TESTING_LEVEL}_ops_level2_epetra "${SRC1}" gTestMain_mpi 3)

 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_epetra_level3.cc)
 add_utest_mpi(${TESTING_LEVEL}_ops_level3_epetra "${SRC1}" gTestMain_mpi 3)
endif()
