
# to include google test fixtures
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/fixtures)

add_subdirectory(meta)
add_subdirectory(serial)

if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_subdirectory(kokkos)
endif()

if(PRESSIO_ENABLE_TPL_MPI)
  add_subdirectory(mpi)
endif()
