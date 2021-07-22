
# to include google test fixtures
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/fixtures)

if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_eigen_vector.cc)
  add_serial_utest(ops_eigen_vector "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_eigen_matrix.cc)
  add_serial_utest(ops_eigen_matrix "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_eigen_diag.cc)
  add_serial_utest(ops_eigen_diag "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_eigen_span.cc)
  add_serial_utest(ops_eigen_span "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_eigen_subspan.cc)
  add_serial_utest(ops_eigen_subspan "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_eigen_level2.cc)
  add_serial_utest(ops_eigen_level2 "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_eigen_level3.cc)
  add_serial_utest(ops_eigen_level3 "${SRC1}")
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_kokkos_vector.cc)
  add_serial_utest_kokkos(ops_kokkos_vector "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_kokkos_matrix.cc)
  add_serial_utest_kokkos(ops_kokkos_matrix "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_kokkos_diag.cc)
  add_serial_utest_kokkos(ops_kokkos_diag "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_kokkos_span.cc)
  add_serial_utest_kokkos(ops_kokkos_span "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_kokkos_subspan.cc)
  add_serial_utest_kokkos(ops_kokkos_subspan "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_kokkos_level2.cc)
  add_serial_utest_kokkos(ops_kokkos_level2 "${SRC1}")

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_kokkos_level3.cc)
  add_serial_utest_kokkos(ops_kokkos_level3 "${SRC1}")
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_tpetra_vector.cc)
 add_utest_mpi(ops_tpetra_vector "${SRC1}" gTestMain_tpetra 3)

 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_tpetra_multi_vector.cc)
 add_utest_mpi(ops_tpetra_multi_vector "${SRC1}" gTestMain_tpetra 3)

 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_tpetra_level2.cc)
 add_utest_mpi(ops_tpetra_level2 "${SRC1}" gTestMain_tpetra 3)

 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_tpetra_level3.cc)
 add_utest_mpi(ops_tpetra_level3 "${SRC1}" gTestMain_tpetra 3)

 #set(SRC2 ${CMAKE_CURRENT_SOURCE_DIR}/ops_teuchos.cc)
 #add_serial_utest(ops_serial_2 "${SRC2}")
endif()
