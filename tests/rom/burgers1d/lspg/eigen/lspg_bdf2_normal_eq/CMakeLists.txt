
configure_file(../../../basis/svd_basis_ncell20_t010_dt001_implicit_bdf2.txt basis.txt COPYONLY)
add_serial_exe_and_test(burgers1d_eigen_lspg_bdf2_gn_normal_eqns rom main.cc "PASSED")
