
set(SRC1
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_wrappers_scalar_compatibility.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_meta_basic.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_meta_detect_typedefs.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/containers_arbitrary_wrappers_static_check.cc
  )
add_serial_utest(containers_meta "${SRC1}")
