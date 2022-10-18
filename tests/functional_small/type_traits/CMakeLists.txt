
set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/detect_typedefs.cc)
add_serial_utest(${TESTING_LEVEL}_type_traits_nested_typedefs ${SRC1})

set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/detect_methods.cc)
add_serial_utest(${TESTING_LEVEL}_type_traits_nested_methods ${SRC1})

#if defined(PRESSIO_ENABLE_TPL_EIGEN) // use Eigen types for testing

if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/all_have_traits.cc)
  set(EXE ${TESTING_LEVEL}_type_traits_all_have_traits_eigen)
  add_serial_utest(${EXE} ${SRC1})
  target_compile_definitions(${EXE} PUBLIC USE_EIGEN)

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/all_have_traits_and_same_scalar.cc)
  set(EXE ${TESTING_LEVEL}_type_traits_all_have_traits_and_same_scalar_eigen)
  add_serial_utest(${EXE} ${SRC1})
  target_compile_definitions(${EXE} PUBLIC USE_EIGEN)

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/traits_eigen.cc)
  add_serial_utest(${TESTING_LEVEL}_type_traits_eigen ${SRC1})
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/all_have_traits.cc)
  set(EXE ${TESTING_LEVEL}_type_traits_all_have_traits_kokkos)
  add_serial_utest(${EXE} ${SRC1})
  target_compile_definitions(${EXE} PUBLIC USE_KOKKOS)

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/all_have_traits_and_same_scalar.cc)
  set(EXE ${TESTING_LEVEL}_type_traits_all_have_traits_and_same_scalar_kokkos)
  add_serial_utest(${EXE} ${SRC1})
  target_compile_definitions(${EXE} PUBLIC USE_KOKKOS)

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

  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/traits_tpetra_block.cc)
  add_serial_utest(${TESTING_LEVEL}_type_traits_tpetra_block ${SRC})
endif()
