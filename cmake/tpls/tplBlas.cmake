
if(PRESSIO_ENABLE_TPL_BLAS)
  message(">> Enabling BLAS since PRESSIO_ENABLE_TPL_BLAS=ON ==> enabling also LAPACK")
  set(PRESSIO_ENABLE_TPL_LAPACK ON)
endif()


if(PRESSIO_ENABLE_TPL_BLAS OR PRESSIO_ENABLE_TPL_TRILINOS)

  if(PRESSIO_ENABLE_UNIT_TESTS OR PRESSIO_ENABLE_TESTS)
    # check if BLAS_ROOT is specified
    if (NOT ${BLAS_ROOT})
      message("")
      message(FATAL_ERROR "BLAS_ROOT not speificed, terminating")
      message("Make sure you set the BLAS_ROOT env var")
    endif()

    cmake_policy(SET CMP0074 NEW)
    find_package(BLAS REQUIRED)
    link_libraries(${BLAS_LIBRARIES})
    message("")
  endif()

endif()
