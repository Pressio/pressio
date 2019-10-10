
option(TPL_ENABLE_MPI "Enable MPI" OFF)
if(TPL_ENABLE_MPI)
  set(HAVE_MPI ON)
  find_package(MPI REQUIRED)

  include_directories(${MPI_CXX_INCLUDE_PATH})
endif()
