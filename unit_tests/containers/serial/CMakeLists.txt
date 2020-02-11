
set(SRC1
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_meta_basic.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_meta_detect_typedefs.cc
  #
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_matrix_dense_serial_eigen_subspan.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_serial_eigen_span.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_multi_vector_eigen_product_eigen_subview_vector_from_mv.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_multi_vector_eigen_view_col_vector.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_multi_vector_eigen_product_eigen_vector_span.cc
  #  
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_matrix_vector_product_eigen.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_matrix_matrix_product_eigen.cc
  #
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_multi_vector_eigen_product_eigen_vector.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_multi_vector_eigen_dot_eigen_vector.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_multi_vector_eigen_dot_multi_vector_eigen.cc
  #
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_serial_printer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_multi_vector_serial_eigen_dynamic_class.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_serial_eigen_dynamic_class.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_matrix_dense_serial_eigen_class.cc
  )
add_serial_utest(containers_serial "${SRC1}")


if(PRESSIO_ENABLE_TPL_TRILINOS)
  set(SRC2
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_serial_teuchos_dense_class.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_matrix_dense_serial_teuchos_class.cc
    )
  add_utest_mpi(containers_teuchos "${SRC2}" gTestMain_mpi 1)
endif()
