
if (PRESSIO_ENABLE_TPL_EIGEN)
  # GN test uses a weighting operator that accepts wrappers
  add_serial_exe_and_test   (gn_neq_dummy_weighting_accept_wrappers_correctness solvers main.cc "PASSED")
  target_compile_definitions(solvers_gn_neq_dummy_weighting_accept_wrappers_correctness PRIVATE USE_GN_NEQ USE_WRAPPERS)
  # GN test uses a weighting operator that accepts native
  add_serial_exe_and_test   (gn_neq_dummy_weighting_accept_native_correctness solvers main.cc "PASSED")
  target_compile_definitions(solvers_gn_neq_dummy_weighting_accept_native_correctness PRIVATE USE_GN_NEQ USE_NATIVE)

  # LM test uses a weighting operator that accepts wrappers
  add_serial_exe_and_test   (lm_neq_dummy_weighting_accept_wrappers_correctness solvers main.cc "PASSED")
  target_compile_definitions(solvers_lm_neq_dummy_weighting_accept_wrappers_correctness PRIVATE USE_LM_NEQ USE_WRAPPERS)
  # LM test uses a weighting operator that accepts native
  add_serial_exe_and_test   (lm_neq_dummy_weighting_accept_native_correctness solvers main.cc "PASSED")
  target_compile_definitions(solvers_lm_neq_dummy_weighting_accept_native_correctness PRIVATE USE_LM_NEQ USE_NATIVE)
endif()
