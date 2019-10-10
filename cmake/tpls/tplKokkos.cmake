
option(TPL_ENABLE_KOKKOS "Enable Kokkos TPL" OFF)
if(TPL_ENABLE_KOKKOS OR TPL_ENABLE_TRILINOS)
  set(HAVE_KOKKOS ON)
  set(KOKKOS_LIB_NAMES  kokkosalgorithms 
                        kokkoscontainers
                        kokkoscore
                        kokkoskernels
)
endif()
