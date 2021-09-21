
configure_file(../../../basis/basis.txt basis.txt COPYONLY)
add_serial_exe_and_test(swe2d_eigen_lspg_crank_nicolson_gn_normal_eqns rom main.cc "PASSED")
