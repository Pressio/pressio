
add_subdirectory(eigen)

if(HAVE_TRILINOS)
  add_subdirectory(epetra)
  add_subdirectory(tpetra)
endif()

if(HAVE_KOKKOS)
  add_subdirectory(kokkos)
endif()

