
configure_file(../../basis/svd_basis_ncell20_t010_dt001_implicit_euler.txt basis.txt COPYONLY)
add_serial_exe_and_test(burgers1d_wl_rom rom main.cc "PASSED")
add_mpi_exe_and_test(burgers1d_wl_rom_mpi rom main.cc 2 "PASSED")
