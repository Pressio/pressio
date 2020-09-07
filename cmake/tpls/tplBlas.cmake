
if(PRESSIO_ENABLE_TPL_BLAS)
  message("> Enabling BLAS since PRESSIO_ENABLE_TPL_BLAS=ON ==> enabling also LAPACK")
  set(PRESSIO_ENABLE_TPL_LAPACK ON)
endif()


if(PRESSIO_ENABLE_TPL_BLAS OR PRESSIO_ENABLE_TPL_TRILINOS)

  if(PRESSIO_ENABLE_UNIT_TESTS OR PRESSIO_ENABLE_TESTS)
    # check if BLAS_DIR is specified
    if (BLAS_DIR)
      message("")
      message("I found BLAS_DIR=${BLAS_DIR}")
      message("If this is not right, reconfigure with: -DBLAS_DIR=<path-to-your-blas>")
    endif()

    find_package( BLAS REQUIRED )
    link_libraries(${BLAS_LIBRARIES})
    message("")
  endif()

endif()
