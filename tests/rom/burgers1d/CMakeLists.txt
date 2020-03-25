
add_subdirectory(eigen)
add_subdirectory(arbitrary)

if(PRESSIO_ENABLE_TPL_TRILINOS)
  add_subdirectory(epetra)
  add_subdirectory(tpetra)
  add_subdirectory(tpetra_block)

endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_subdirectory(kokkos)
  # add_subdirectory(wls)
endif()
