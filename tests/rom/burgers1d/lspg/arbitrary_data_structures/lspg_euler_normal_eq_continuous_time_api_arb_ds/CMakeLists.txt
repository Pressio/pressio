
configure_file(../../../basis/svd_basis_ncell20_t010_dt001_implicit_euler.txt basis.txt COPYONLY)
add_serial_exe_and_test(burgers1d_eigen_lspg_euler_gn_normal_eqns_continuous_time_api_arbds rom main.cc "PASSED")
