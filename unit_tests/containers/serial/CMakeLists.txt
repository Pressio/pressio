
set(SRC1
  #
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_serial_eigen_span.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_matrix_dense_serial_eigen_subspan.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_matrix_dense_serial_eigen_diag_expressions.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_matrix_dense_serial_eigen_diag_ops.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_multi_vector_serial_eigen_dynamic_class.cc
  #
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_serial_printer.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_serial_eigen_dynamic_class.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_matrix_dense_serial_eigen_class.cc
  )
add_serial_utest(containers_serial "${SRC1}")

if(PRESSIO_ENABLE_TPL_PYBIND11)
  set(SRC2 ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_pybind.cc)
  add_serial_utest(containers_pybind11 "${SRC2}")
endif()

# if(PRESSIO_ENABLE_TPL_TRILINOS)
#   set(SRC2
#     ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_serial_teuchos_dense_class.cc
#     ${CMAKE_CURRENT_SOURCE_DIR}/containers_matrix_dense_serial_teuchos_class.cc
#     )
#   add_utest_mpi(containers_teuchos "${SRC2}" gTestMain_mpi 1)
# endif()
