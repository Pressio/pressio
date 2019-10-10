

add_subdirectory(serial)

if(HAVE_MPI)
  add_subdirectory(mpi)
endif()
