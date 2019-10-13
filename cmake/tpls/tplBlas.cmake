
option(TPL_ENABLE_BLAS "Enable Blas TPL" OFF)

if(TPL_ENABLE_BLAS OR TPL_ENABLE_TRILINOS)
  set(HAVE_BLAS ON)

  # check if BLAS_DIR is specified
  if (BLAS_DIR)
    message("")
    message(
      "I found BLAS_DIR=${BLAS_DIR}.
If this is not right, then reconfigure with: -DBLAS_DIR=<path-to-your-blas>")
  endif()

  find_package( BLAS REQUIRED )
  link_libraries(${BLAS_LIBRARIES})
  message("")

endif()
