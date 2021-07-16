
set(SRC1
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_meta_detect_typedefs.cc
  if(PRESSIO_ENABLE_TPL_EIGEN)
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_is_are_wrappers.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_wrappers_scalar_compatibility.cc
  endif()
  if(PRESSIO_ENABLE_TPL_TRILINOS)
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_meta_basic.cc
  endif()
  )
add_serial_utest(containers_meta ${SRC1})

add_serial_utest(containers_serial_printer
	${CMAKE_CURRENT_SOURCE_DIR}/containers_serial_printer.cc)

if(PRESSIO_ENABLE_TPL_TRILINOS)
add_utest_mpi(containers_tpetra_maps_import
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_tpetra_maps_import.cc gTestMain_tpetra 3 )
endif()
