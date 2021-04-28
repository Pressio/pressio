

if(PRESSIO_ENABLE_TPL_LAPACK)
  message(">> Enabling LAPACK since PRESSIO_ENABLE_TPL_LAPACK=ON ==> enabling also BLAS")
  set(PRESSIO_ENABLE_TPL_BLAS ON)
endif()

if(PRESSIO_ENABLE_TPL_LAPACK)

  if(PRESSIO_ENABLE_UNIT_TESTS OR PRESSIO_ENABLE_TESTS)
    # check if LAPACK_ROOT is specified
    if (NOT ${LAPACK_ROOT})
      message("")
      message(FATAL_ERROR "LAPACK_ROOT not speificed, terminating")
      message("Make sure you set the LAPACK_ROOT env var")
    endif()

    cmake_policy(SET CMP0074 NEW)
    find_package( LAPACK REQUIRED )
    link_libraries(${LAPACK_LIBRARIES})
    message("LAPLIBS=${LAPACK_LIBRARIES}")
  endif()

endif()
