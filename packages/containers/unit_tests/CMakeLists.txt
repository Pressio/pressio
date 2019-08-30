
#------------------------------------------------------------------------
# to include the lumped headers to have access to containers classes
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../..)
# to include google test fixtures
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/fixtures)

#------------------------------------------------------------------------
# include unit tests subfolders
#------------------------------------------------------------------------
ADD_SUBDIRECTORY(serial)

IF(HAVE_MPI)
  ADD_SUBDIRECTORY(mpi)
ENDIF()

IF(HAVE_KOKKOS)
  ADD_SUBDIRECTORY(kokkos)
ENDIF()
