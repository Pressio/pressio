
option(PRESSIO_ENABLE_TPL_LAPACK "Enable Lapack TPL" OFF)
if(PRESSIO_ENABLE_TPL_LAPACK)

  if(PRESSIO_ENABLE_UNIT_TESTS OR PRESSIO_ENABLE_TESTS)
    # check if LAPACK_DIR is specified
    if (LAPACK_DIR)
      message("")
      message("I found LAPACK_DIR=${LAPACK_DIR}.")
      message("If this is not right, reconfigure with: -DLAPACK_DIR=<path-to-your-lapack>")
    endif()

    find_package( LAPACK REQUIRED )
    link_libraries(${LAPACK_LIBRARIES})
    message("")
  endif()

endif()
