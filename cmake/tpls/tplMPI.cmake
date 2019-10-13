
option(PRESSIO_ENABLE_TPL_MPI "Enable MPI" OFF)
if(PRESSIO_ENABLE_TPL_MPI)
  if(PRESSIO_ENABLE_UNIT_TESTS OR PRESSIO_ENABLE_TESTS)
    find_package(MPI REQUIRED)
    include_directories(${MPI_CXX_INCLUDE_PATH})
  endif()
endif()
