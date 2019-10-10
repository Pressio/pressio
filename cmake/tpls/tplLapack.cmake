
option(TPL_ENABLE_LAPACK "Enable Lapack TPL" OFF)
if(TPL_ENABLE_LAPACK)
  set(HAVE_LAPACK ON)
  find_package( LAPACK REQUIRED )
  link_libraries(${LAPACK_LIBRARIES})
endif()
