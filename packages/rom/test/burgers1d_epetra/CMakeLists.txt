
IF(HAVE_TRILINOS)
ADD_SUBDIRECTORIES(
  lspg_euler_tsqr
  lspg_euler_teuchos_state
  lspg_euler_normal_eq
  lspg_euler_normal_eq_masked
  lspg_euler_normal_eq_reduced_no_mask
  lspg_bdf2_normal_eq
  galerkin_explicit_euler
  lspg_euler_modified_gram_schmidt_qr
  )
ENDIF()
