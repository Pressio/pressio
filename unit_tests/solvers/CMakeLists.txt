
add_subdirectory(meta)
add_subdirectory(serial)

if(PRESSIO_ENABLE_TPL_MPI)
  add_subdirectory(mpi)
endif()
if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_subdirectory(kokkos)
endif()
