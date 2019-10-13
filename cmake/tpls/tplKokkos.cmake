
option(PRESSIO_ENABLE_TPL_KOKKOS "Enable Kokkos TPL" OFF)

if(PRESSIO_ENABLE_TPL_KOKKOS)
  message("Enabling Kokkos")

  if(PRESSIO_ENABLE_UNIT_TESTS OR PRESSIO_ENABLE_TESTS)
    set(KOKKOS_LIB_NAMES
      kokkosalgorithms  kokkoscontainers
      kokkoscore kokkoskernels)
  endif()

endif()
