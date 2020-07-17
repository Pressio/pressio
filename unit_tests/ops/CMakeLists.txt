
# to include google test fixtures
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../containers/fixtures)

add_subdirectory(serial)

if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_subdirectory(kokkos)
endif()

if(PRESSIO_ENABLE_TPL_MPI AND PRESSIO_ENABLE_TPL_TRILINOS)
  add_subdirectory(trilinos)
endif()
