
set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/detect_typedefs.cc)
add_serial_utest(${TESTING_LEVEL}_type_traits_nested_typedefs ${SRC1})

set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/detect_methods.cc)
add_serial_utest(${TESTING_LEVEL}_type_traits_nested_methods ${SRC1})

set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/scalar_compatibility.cc)
add_serial_utest(${TESTING_LEVEL}_type_traits_scalar_compatibility ${SRC1})

if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/traits_eigen.cc)
  add_serial_utest(${TESTING_LEVEL}_type_traits_eigen ${SRC1})
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/traits_kokkos.cc)
  add_serial_utest_kokkos(${TESTING_LEVEL}_type_traits_kokkos ${SRC1})
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/traits_teuchos.cc)
  add_serial_utest(${TESTING_LEVEL}_type_traits_teuchos ${SRC})

  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/traits_epetra.cc)
  add_serial_utest(${TESTING_LEVEL}_type_traits_epetra ${SRC})

  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/traits_tpetra.cc)
  add_serial_utest(${TESTING_LEVEL}_type_traits_tpetra ${SRC})
endif()
