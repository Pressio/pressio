
IF(HAVE_TRILINOS AND HAVE_ANASAZI_TSQR)
ADD_SUBDIRECTORIES(
  lspg_euler_tsqr
  lspg_euler_teuchos_state
  )
ENDIF()

IF(HAVE_TRILINOS)
ADD_SUBDIRECTORIES(
  lspg_euler_normal_eq
  lspg_euler_normal_eq_masked
  lspg_euler_normal_eq_reduced_no_mask
  lspg_bdf2_normal_eq
  galerkin_explicit_euler
  )
ENDIF()
