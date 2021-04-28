
if(PRESSIO_ENABLE_TPL_KOKKOS)
  message(">> Enabling Kokkos since PRESSIO_ENABLE_TPL_KOKKOS=ON")
endif()


if(PRESSIO_ENABLE_TPL_KOKKOS)

  if(PRESSIO_ENABLE_UNIT_TESTS OR PRESSIO_ENABLE_TESTS)

    # # when trilinos is also enabled it links kokkos too, see tplTrilinos.cmake
    if(NOT PRESSIO_ENABLE_TPL_TRILINOS)
      # if kokkos is used as standalone lib, then we are more specific
      # user needs to defined: KOKKOS_ROOT_DIR and KOKKOS_KERNELS_ROOT_DIR
      if (NOT KOKKOS_ROOT OR NOT KOKKOS_KERNELS_ROOT)
	message(FATAL_ERROR "Missing KOKKOS_ROOT. KOKKOS needs:
          -D KOKKOS_ROOT=<full-path-to-kokkos-installation>
          -D KOKKOS_KERNELS_ROOT=<full-path-to-kokkos-kernels-installation>
          ")
      endif()

      set(KOKKOS_LIB_NAMES kokkoscontainers kokkoscore kokkoskernels)

      include_directories(${KOKKOS_ROOT}/include ${KOKKOS_KERNELS_ROOT}/include)

      link_directories(${KOKKOS_ROOT}/lib ${KOKKOS_ROOT}/lib64
	${KOKKOS_KERNELS_ROOT}/lib ${KOKKOS_KERNELS_ROOT}/lib64)

      link_libraries(${KOKKOS_LIB_NAMES})
    endif()

  endif()

endif()
