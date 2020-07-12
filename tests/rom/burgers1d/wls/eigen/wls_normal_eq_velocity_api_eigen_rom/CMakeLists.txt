
configure_file(../../../basis/svd_basis_ncell20_t010_dt001_implicit_euler.txt basis_euler.txt COPYONLY)
configure_file(../../../basis/svd_basis_ncell20_t010_dt001_implicit_bdf2.txt basis_bdf2.txt COPYONLY)

add_serial_exe_and_test(burgers1d_eigen_wls_euler_gn_normal_eqns_eigen rom main.cc "PASSED")
