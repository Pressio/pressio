
# to include google test fixtures
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/fixtures)

if(PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_utest(
    ${TESTING_LEVEL}_ops_eigen
    ops_eigen_vector.cc
    ops_eigen_matrix.cc
    ops_eigen_diag.cc
    ops_eigen_column.cc
    ops_eigen_span.cc
    ops_eigen_subspan.cc
    ops_eigen_level2.cc
    ops_eigen_level3.cc
  )
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_serial_utest_kokkos(${TESTING_LEVEL}_ops_vector_kokkos
    ops_kokkos_vector.cc
    ops_kokkos_matrix.cc
    ops_kokkos_diag.cc
    ops_kokkos_span.cc
    ops_kokkos_subspan.cc
    ops_kokkos_level2.cc
    ops_kokkos_level3.cc
  )
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

 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_tpetra_column_expr.cc)
 add_utest_mpi(${TESTING_LEVEL}_ops_column_expr_tpetra "${SRC1}" gTestMain_tpetra 3)


 # TPETRA BLOCK
 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_tpetra_block_vector.cc)
 add_utest_mpi(${TESTING_LEVEL}_ops_vector_tpetra_block "${SRC1}" gTestMain_tpetra 3)

 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_tpetra_block_column_expr.cc)
 add_utest_mpi(${TESTING_LEVEL}_ops_column_expr_tpetra_block "${SRC1}" gTestMain_tpetra 3)

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

 # TEUCHOS
 set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ops_teuchos_vector.cc)
 add_serial_utest(${TESTING_LEVEL}_ops_vector_teuchos "${SRC1}")
endif()
