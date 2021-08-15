
configure_file(../../../basis/svd_basis_ncell20_t010_dt001_implicit_euler.txt basis.txt COPYONLY)
add_mpi_exe_and_test(burgers1d_tpetra_block_lspg_kokkos_rom_euler_gn_qr rom main.cc 1 "PASSED")
