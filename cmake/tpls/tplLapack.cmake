
option(TPL_ENABLE_LAPACK "Enable Lapack TPL" OFF)
if(TPL_ENABLE_LAPACK)
  set(HAVE_LAPACK ON)

  # check if LAPACK_DIR is specified
  if (LAPACK_DIR)
    message("")
    message(
      "I found LAPACK_DIR=${LAPACK_DIR}.
If this is not right, then reconfigure with: -DLAPACK_DIR=<path-to-your-lapack>")
  endif()

  find_package( LAPACK REQUIRED )
  link_libraries(${LAPACK_LIBRARIES})
  message("")

endif()
