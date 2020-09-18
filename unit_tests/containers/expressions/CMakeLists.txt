
set(exeRootName containers_expressions)

set(SRC1
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_serial_eigen_span.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_matrix_dense_serial_eigen_subspan.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_matrix_dense_serial_eigen_diag_expressions.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_expressions_are_static_eigen.cc)
add_serial_utest(${exeRootName}_serial "${SRC1}")


if(PRESSIO_ENABLE_TPL_KOKKOS)
add_serial_utest_kokkos(${exeRootName}_dense_kokkos_diag
  	${CMAKE_CURRENT_SOURCE_DIR}/containers_matrix_dense_kokkos_diag_expression.cc)

add_serial_utest_kokkos(${exeRootName}_are_static_kokkos
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_expressions_are_static_kokkos.cc)
endif()
