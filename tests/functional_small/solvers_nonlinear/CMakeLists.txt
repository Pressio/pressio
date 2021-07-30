
set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/concepts.cc)
add_serial_utest(${TESTING_LEVEL}_solvers_nonlinear_concepts ${SRC1})

if (PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/tolerances.cc)
  add_serial_utest(${TESTING_LEVEL}_solvers_nonlinear_tolerances ${SRC1})

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/operators.cc)
  add_serial_utest(${TESTING_LEVEL}_solvers_nonlinear_operators ${SRC1})

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/correctors.cc)
  add_serial_utest(${TESTING_LEVEL}_solvers_nonlinear_correctors ${SRC1})

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/ls_weighted_gn_lm_gt.cc)
  set(TESTNAME ${TESTING_LEVEL}_solvers_nonlinear_weighted_least_squares_gn)
  add_serial_utest(${TESTNAME} ${SRC1})

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/iteratively_reweighted_gn.cc)
  set(TESTNAME ${TESTING_LEVEL}_solvers_nonlinear_irwls_gn)
  add_serial_utest(${TESTNAME} ${SRC1})
endif()
