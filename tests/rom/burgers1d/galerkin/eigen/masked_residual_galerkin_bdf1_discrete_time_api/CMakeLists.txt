
configure_file(../../../basis/svd_basis_ncell20_t010_dt001_implicit_euler.txt basis.txt COPYONLY)
set(myname burgers1d_masked_residual_galerkin_bdf1_discrete_time_api)
add_serial_exe_and_test(${myname} rom main.cc "PASSED")
