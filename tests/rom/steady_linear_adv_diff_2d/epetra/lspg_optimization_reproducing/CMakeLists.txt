
configure_file(../../basis/basis_row189_col5.txt basis.txt COPYONLY)
add_mpi_exe_and_test(steady_lin_adv_diff2d_epetra_lspg_opt_reproducing rom main.cc 1 "PASSED")
