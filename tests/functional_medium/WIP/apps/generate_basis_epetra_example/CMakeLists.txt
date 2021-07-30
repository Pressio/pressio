
configure_file(gold_basis.txt gold_basis.txt COPYONLY)
add_mpi_exe_and_test(steadyAdvDiff2d_fom_gen_basis_epetra_example apps main.cc 1 "PASSED")
