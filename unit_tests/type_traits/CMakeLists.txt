
set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/detect_typedefs.cc)
add_serial_utest(type_traits_nested_typedefs ${SRC1})

if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/traits_eigen.cc)
  add_serial_utest(type_traits_eigen ${SRC1})
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_serial_utest_kokkos(type_traits_kokkos ${CMAKE_CURRENT_SOURCE_DIR}/traits_kokkos.cc )
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/traits_teuchos.cc)
  add_serial_utest(type_traits_teuchos ${SRC})

  ### EPETRA ###
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/traits_epetra.cc)
  add_serial_utest(type_traits_epetra ${SRC})

  ### TPETRA ###
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/traits_tpetra.cc)
  add_serial_utest(type_traits_tpetra ${SRC})

  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/traits_tpetra_block.cc)
  add_serial_utest(type_traits_tpetra_block ${SRC})
endif()
