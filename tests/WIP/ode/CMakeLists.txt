
add_subdirectory(burgers1d)
add_subdirectory(swe2d)

if(PRESSIO_ENABLE_TPL_EIGEN)
  add_subdirectory(advance_n_steps_arbitrary_dt_setter_with_collector)
  add_subdirectory(integrate_to_time_with_time_step_recovery)
  add_subdirectory(crank_nicolson_correctness)
endif()
