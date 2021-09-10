
if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/eigen_traits.cc)
  add_serial_utest(${TESTING_LEVEL}_expressions_eigen_traits ${SRC})

  set(SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/eigen_span.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/eigen_subspan.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/eigen_diag.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/eigen_asdiagonalmatrix.cc
    )
  add_serial_utest(${TESTING_LEVEL}_expressions_eigen_numerics ${SRC})
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/kokkos_traits.cc)
  add_serial_utest_kokkos(${TESTING_LEVEL}_expressions_kokkos_traits ${SRC})

  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/kokkos_span.cc)
  add_serial_utest_kokkos(${TESTING_LEVEL}_expressions_kokkos_span ${SRC})

  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/kokkos_subspan.cc)
  add_serial_utest_kokkos(${TESTING_LEVEL}_expressions_kokkos_subspan ${SRC})

  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/kokkos_diag.cc)
  add_serial_utest_kokkos(${TESTING_LEVEL}_expressions_kokkos_diag ${SRC})

  add_serial_utest_kokkos(${TESTING_LEVEL}_kokkos_strided_subview 
    ${CMAKE_CURRENT_SOURCE_DIR}/kokkos_strided_subview.cc )
endif()
