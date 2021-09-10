
configure_file(../../../basis/svd_basis_ncell20_t010_dt001_implicit_euler.txt basis.txt COPYONLY)
add_mpi_exe_and_test(
  burgers1d_epetra_masked_velocity_galerkin_explicit_euler_continuous_time_api
  rom main.cc 3 "PASSED")
