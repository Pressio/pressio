
if (PRESSIO_ENABLE_TPL_EIGEN)
  add_serial_exe_and_test   (gn_neq_trivial_weighting_correctness solvers_nonlinear main.cc "PASSED")
  target_compile_definitions(solvers_nonlinear_gn_neq_trivial_weighting_correctness PRIVATE USE_GN_NEQ)

  add_serial_exe_and_test   (lm_neq_trivial_weighting_correctness solvers_nonlinear main.cc "PASSED")
  target_compile_definitions(solvers_nonlinear_lm_neq_trivial_weighting_correctness PRIVATE USE_LM_NEQ)
endif()
