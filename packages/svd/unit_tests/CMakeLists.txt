
#------------------------------------------------------------------------
# to include the lumped headers to have access to core classes
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../..)

#------------------------------------------------------------------------
# include unit tests subfolders
#------------------------------------------------------------------------
ADD_SUBDIRECTORY(serial)

IF(HAVE_MPI)
  ADD_SUBDIRECTORY(mpi)
ENDIF()