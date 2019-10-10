
option(TPL_ENABLE_BLAS "Enable Blas TPL" OFF)

if(TPL_ENABLE_BLAS OR TPL_ENABLE_TRILINOS)
  set(HAVE_BLAS ON)
  find_package( BLAS REQUIRED )
  link_libraries(${BLAS_LIBRARIES})
endif()
