
#------------------------------------------------------------------------
# to include the lumped headers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../..)

IF(HAVE_TRILINOS)
ADD_SUBDIRECTORY(epetra)
ADD_SUBDIRECTORY(generate_basis_epetra_example)
ENDIF()
