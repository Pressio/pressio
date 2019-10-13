
add_subdirectory(eigen)

if(PRESSIO_ENABLE_TPL_TRILINOS)
  add_subdirectory(epetra)
  add_subdirectory(tpetra)
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_subdirectory(kokkos)
endif()
