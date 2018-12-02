

#to include the lumped headers
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../..)


IF(HAVE_TRILINOS)
ADD_SUBDIRECTORIES(
	# svd_epetra_wrapper
	#svd_trilinos_anasazi
)
ENDIF()
