
configure_file(U20_ncell50_t35_dt001.txt basis.txt COPYONLY)
add_mpi_exe_and_test(burgers1d_epetra_galerkin_explicit_euler rom main.cc 2 "PASSED")
