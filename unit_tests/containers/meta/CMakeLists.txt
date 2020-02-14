
set(SRC1
  ${CMAKE_CURRENT_SOURCE_DIR}/wrappers_scalar_compatibility.cc
  )
add_serial_utest(containers_meta "${SRC1}")
