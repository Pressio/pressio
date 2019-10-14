
if(PRESSIO_ENABLE_TPL_TRILINOS)
  add_utest_mpi(rom_default_lspg_epetra rom_default_lspg_epetra.cc gTestMain_mpi 1)
  add_utest_mpi(rom_states_fom_epetra rom_states_fom_epetra.cc gTestMain_mpi 2)
endif()
