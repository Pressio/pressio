

add_subdirectory(arbitrary_data_structures)
add_subdirectory(eigen)

if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_subdirectory(kokkos)
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
  add_subdirectory(epetra)
  add_subdirectory(tpetra)
  add_subdirectory(tpetra_block)
endif()
