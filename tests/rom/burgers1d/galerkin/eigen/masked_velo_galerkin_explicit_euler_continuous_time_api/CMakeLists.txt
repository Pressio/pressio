
configure_file(../../../basis/svd_basis_ncell20_t010_dt001_implicit_euler.txt basis.txt COPYONLY)
add_serial_exe_and_test(
  burgers1d_masked_velo_galerkin_explicit_euler_continuous_time_api
  rom main.cc "PASSED")
