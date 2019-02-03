
#------------------------------------------------------------------------
# to include the lumped headers to have access to core classes
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/../..)


IF(HAVE_TRILINOS AND HAVE_ANASAZI_TSQR)
ADD_SUBDIRECTORIES(
 burgers1d_epetra_lspg_euler_tsqr
 burgers1d_epetra_lspg_euler_teuchos_state
)
ENDIF()

ADD_SUBDIRECTORIES(
  burgers1d_epetra_lspg_euler
  burgers1d_epetra_lspg_bdf2
)
