
if(PRESSIO_ENABLE_TPL_EIGEN)
  add_subdirectory(eigen)
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
  add_subdirectory(tpetra)
  add_subdirectory(tpetra_block)
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_subdirectory(kokkos)
endif()
