

if(PRESSIO_ENABLE_TPL_TRILINOS)
  add_subdirectory(tpetra)
  add_subdirectory(tpetra_hyperreduction)

endif()

#if(PRESSIO_ENABLE_TPL_KOKKOS)
#  add_subdirectory(kokkos)
#endif()
