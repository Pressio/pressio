
if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/solvers_linear_eigen.cc)
  add_serial_utest(${TESTING_LEVEL}_solvers_linear_eigen ${SRC1})
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/solvers_linear_kokkos.cc)
  add_serial_utest_kokkos(${TESTING_LEVEL}_solvers_linear_kokkos ${SRC1})
endif()
