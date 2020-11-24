
if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC1
    ${CMAKE_CURRENT_SOURCE_DIR}/ops_multi_vector_eigen_product_eigen_vector_span.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ops_multi_vector_eigen_product_eigen_vector.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ops_multi_vector_eigen_dot_eigen_vector.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ops_multi_vector_eigen_dot_multi_vector_eigen.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ops_matrix_vector_product_eigen.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ops_matrix_matrix_product_eigen.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ops_vector_serial_eigen_dynamic_class.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ops_matrix_dense_serial_eigen_class.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ops_matrix_dense_serial_eigen_diag.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ops_asdiagonalmatrix_eigen_product_eigen_matrix.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ops_vector_eigen_pow.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/ops_vector_eigen_abs.cc
    )
  add_serial_utest(ops_serial "${SRC1}")
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
set(SRC2
  ${CMAKE_CURRENT_SOURCE_DIR}/ops_teuchos_matrix_eigen_vector_product.cc
  )
add_serial_utest(ops_serial_2 "${SRC2}")
endif()
