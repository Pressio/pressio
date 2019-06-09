
TRIBITS_ADD_EXECUTABLE(
ks1d_rom_eigen_lspg_bdf2_gn_qr
  SOURCES main.cc
  )

configure_file(basis/svd_basis_ncell20_t010_dt001_implicit_bdf2.txt basis.txt COPYONLY)
