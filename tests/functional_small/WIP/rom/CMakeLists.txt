
add_subdirectory(serial)

if(PRESSIO_ENABLE_TPL_MPI)
  add_subdirectory(mpi)
endif()
