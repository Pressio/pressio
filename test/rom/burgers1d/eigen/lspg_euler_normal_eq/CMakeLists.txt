
# TRIBITS_ADD_EXECUTABLE_AND_TEST(
#   burgers1d_rom_eigen_lspg_euler_gn_normal_eqns
#   SOURCES main.cc
#   NUM_MPI_PROCS 1
#   PASS_REGULAR_EXPRESSION "PASSED"
#   )

set(exeName burgers1d_rom_eigen_lspg_euler_gn_normal_eqns)
add_executable(${exeName} main.cc)
add_test(NAME ${exeName} COMMAND ${exeName})

configure_file(../../basis/svd_basis_ncell20_t010_dt001_implicit_euler.txt basis.txt COPYONLY)
