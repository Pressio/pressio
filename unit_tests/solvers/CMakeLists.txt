
add_subdirectory(meta)
add_subdirectory(serial)

if(HAVE_MPI)
  add_subdirectory(mpi)
endif()
if(HAVE_KOKKOS)
  add_subdirectory(kokkos)
endif()
