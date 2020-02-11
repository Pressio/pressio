
if(PRESSIO_ENABLE_TPL_KOKKOS AND PRESSIO_ENABLE_TPL_TRILINOS)
  add_subdirectory(lspg_euler_normal_eq)
  # add_subdirectory(wls_euler_normal_eq)
  # add_subdirectory(wls_bdf2_normal_eq)
endif()
