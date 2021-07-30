
configure_file(U20_ncell50_t35_dt001.txt basis.txt COPYONLY)
add_serial_exe_and_test(burgers1d_default_galerkin_velocity_api_arbds rom main.cc "PASSED")
