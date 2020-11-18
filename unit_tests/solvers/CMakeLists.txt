
if (PRESSIO_ENABLE_TPL_EIGEN)
  add_subdirectory(api_checks)

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/solvers_nonlinear_tolerances.cc)
  add_serial_utest(solvers_nonlinear_tolerances ${SRC1})
endif()
