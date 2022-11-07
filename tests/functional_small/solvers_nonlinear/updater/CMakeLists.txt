
if (PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/basic_backtrack.cc)
  add_serial_utest(${TESTING_LEVEL}_basic_backtrack ${SRC1})
endif()
