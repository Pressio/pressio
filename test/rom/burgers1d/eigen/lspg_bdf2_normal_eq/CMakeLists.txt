
set(exeName burgers1d_rom_eigen_lspg_bdf2_gn_normal_eqns)
add_executable(${exeName} main.cc)
add_test(NAME ${exeName} COMMAND ${exeName})

configure_file(../../basis/svd_basis_ncell20_t010_dt001_implicit_bdf2.txt basis.txt COPYONLY)
