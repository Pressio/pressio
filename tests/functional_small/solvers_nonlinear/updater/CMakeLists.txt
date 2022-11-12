
if (PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/backtrack_strictly_decreasing_objective)
  add_serial_utest(${TESTING_LEVEL}_backtrack_strictly_decreasing_objective  ${SRC1})
endif()
