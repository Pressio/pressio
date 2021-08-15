
configure_file(../../../basis/svd_basis_ncell20_t010_dt001_implicit_euler.txt basis.txt COPYONLY)
add_serial_exe_and_test(burgers1d_rom_kokkos_lspg_euler_gn_normal_eqns rom main.cc "PASSED")
