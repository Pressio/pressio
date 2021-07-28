
set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/concepts.cc)
add_serial_utest(solvers_nonlinear_concepts ${SRC1})

if (PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/tolerances.cc)
  add_serial_utest(solvers_nonlinear_tolerances ${SRC1})
endif()
