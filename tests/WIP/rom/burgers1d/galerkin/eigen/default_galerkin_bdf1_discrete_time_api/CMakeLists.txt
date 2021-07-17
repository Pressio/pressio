
configure_file(../../../basis/svd_basis_ncell20_t010_dt001_implicit_euler.txt basis.txt COPYONLY)

# test the case where the galerkin jacobian is specified
set(myname burgers1d_default_galerkin_bdf1_discrete_time_api_with_specified_gal_jacobian)
add_serial_exe_and_test(${myname} rom main.cc "PASSED")
target_compile_definitions(rom_${myname} PRIVATE PASS_GALERKIN_JACOBIAN)

# test the case where the galerkin jacobian is chosen by pressio
set(myname burgers1d_default_galerkin_bdf1_discrete_time_api_with_default_gal_jacobian)
add_serial_exe_and_test(${myname} rom main.cc "PASSED")
