

set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/containers_matrix_dense_serial_eigen_class.cc)
add_serial_utest(containers_matrix_serial "${SRC1}")


if(PRESSIO_ENABLE_TPL_TRILINOS)
  add_serial_utest(containers_matrix_sparse_tpetra_traits
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_sparse_matrix_tpetra_traits.cc)
endif()

