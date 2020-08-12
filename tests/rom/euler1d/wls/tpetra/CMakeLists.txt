
configure_file(../../basis/euler_basis_ncell_500_dt_5em4_jrk4_k61.txt basis.txt COPYONLY)
add_mpi_exe_and_test(euler1d_tpetra apps main.cc 2 "PASSED")
