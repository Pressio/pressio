

set(SRC1
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers_admissible_system_for_nonlinear_solve.cc
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers_nonadmissible_system_for_nonlinear_solve.cc
  )
add_serial_utest(solvers_meta_serial "${SRC1}")

set(SRC2
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers_gn_normal_eq.cc
  )
add_serial_utest(solvers_meta_serial_gn_neq "${SRC2}")

set(SRC3
  ${CMAKE_CURRENT_SOURCE_DIR}/solvers_gn_with_qr.cc
  )
add_serial_utest(solvers_meta_serial_gn_qr "${SRC3}")
