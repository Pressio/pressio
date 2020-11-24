
set(exeRootName containers_expressions)

if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC1
    ${CMAKE_CURRENT_SOURCE_DIR}/eigen_span.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/eigen_subspan.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/eigen_diag.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/eigen_asdiagonalmatrix.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/expressions_are_static_eigen.cc
    )
  add_serial_utest(${exeRootName}_serial "${SRC1}")
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
add_serial_utest_kokkos(${exeRootName}_are_static_kokkos
  ${CMAKE_CURRENT_SOURCE_DIR}/expressions_are_static_kokkos.cc)

add_serial_utest_kokkos(${exeRootName}_kokkos_span
  ${CMAKE_CURRENT_SOURCE_DIR}/kokkos_span.cc)

add_serial_utest_kokkos(${exeRootName}_kokkos_subspan
  ${CMAKE_CURRENT_SOURCE_DIR}/kokkos_subspan.cc)

add_serial_utest_kokkos(${exeRootName}_kokkos_diag
  ${CMAKE_CURRENT_SOURCE_DIR}/kokkos_diag.cc)
endif()
