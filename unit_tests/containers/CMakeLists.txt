
# to include google test fixtures
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/fixtures)

#add_subdirectory(serial)

if(HAVE_MPI)
  add_subdirectory(mpi)
endif()

# if(HAVE_KOKKOS)
#   add_subdirectory(kokkos)
# endif()
