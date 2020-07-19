
option(PRESSIO_ENABLE_TPL_KOKKOS "Enable Kokkos TPL" OFF)

if(PRESSIO_ENABLE_TPL_KOKKOS)
  message("Enabling Kokkos")

  if(PRESSIO_ENABLE_UNIT_TESTS OR PRESSIO_ENABLE_TESTS)

    if(PRESSIO_ENABLE_TPL_TRILINOS)
      # when trilinos is also enabled, look for Kokkos inside there
      # so use the libraries that are built by trilinos
      # and don't worry about include and lib_dir since these are set
      # by the tplTrilinos.cmake
      set(KOKKOS_LIB_NAMES kokkosalgorithms  kokkoscontainers kokkoscore kokkoskernels)

    else()
      # if kokkos is used as standalone lib, then we are more specific
      # user needs to defined: KOKKOS_ROOT_DIR and KOKKOS_KERNELS_ROOT_DIR
      if (NOT KOKKOS_ROOT or KOKKOS_KERNELS_ROOT)
	message(FATAL_ERROR
	  "Missing KOKKOS_ROOT. KOKKOS needs:
          -D KOKKOS_ROOT=<full-path-to-kokkos-installation>
          -D KOKKOS_KERNELS_ROOT=<full-path-to-kokkos-kernels-installation>
          ")
      endif()

      set(KOKKOS_LIB_NAMES kokkoscontainers kokkoscore kokkoskernels)

      include_directories(${KOKKOS_ROOT}/include ${KOKKOS_KERNELS_ROOT}/include)

      include_directories(${KOKKOS_ROOT}/lib ${KOKKOS_ROOT}/lib64
	${KOKKOS_KERNELS_ROOT}/lib ${KOKKOS_KERNELS_ROOT}/lib64)

      link_libraries(${KOKKOS_LIB_NAMES})
    endif()

  endif()

endif()
