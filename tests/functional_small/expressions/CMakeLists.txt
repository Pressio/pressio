# to include google test fixtures
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/fixtures)

if(PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_utest(
    ${TESTING_LEVEL}_expressions_eigen
    eigen_span.cc
    eigen_subspan.cc
    eigen_diag.cc
    eigen_column.cc
    )
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_serial_utest_kokkos(
    ${TESTING_LEVEL}_expressions_kokkos
    kokkos_span.cc
    kokkos_subspan.cc
    kokkos_diag.cc
    kokkos_strided_subview.cc
  )
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/tpetra_column.cc)
  add_utest_mpi(${TESTING_LEVEL}_expressions_tpetra_column gTestMain_tpetra 3 "${SRC}")

  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/tpetra_block_column.cc)
  add_utest_mpi(${TESTING_LEVEL}_expressions_tpetra_block_column gTestMain_tpetra 3 "${SRC}")
endif()
