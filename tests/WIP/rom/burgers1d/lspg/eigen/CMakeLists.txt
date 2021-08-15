
add_subdirectory(lspg_euler_normal_eq)
add_subdirectory(lspg_crank_nicolson_normal_eq)
add_subdirectory(lspg_bdf2_normal_eq)
add_subdirectory(lspg_euler_normal_eq_user_dt_setter)
add_subdirectory(lspg_euler_normal_eq_discrete_time_api)
add_subdirectory(lspg_euler_normal_eq_discrete_time_api_user_defined_mapping)

if(PRESSIO_ENABLE_TPL_TORCH)
  if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    add_subdirectory(lspg_euler_normal_eq_discrete_time_api_libtorch_mapping)
  endif()
endif()
