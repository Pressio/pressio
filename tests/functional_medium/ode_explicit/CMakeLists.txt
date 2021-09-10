
set(ROOTNAME ${TESTING_LEVEL}_ode_explicit)

if(PRESSIO_ENABLE_TPL_EIGEN)
add_serial_exe_and_test(burgers1d_euler_eigen ${ROOTNAME} burgers1d_euler_eigen.cc "PASSED")
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
add_serial_exe_and_test(burgers1d_euler_kokkos ${ROOTNAME} burgers1d_euler_kokkos.cc "PASSED")
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
add_mpi_exe_and_test(burgers1d_euler_epetra ${ROOTNAME} burgers1d_euler_epetra.cc 4 "PASSED")
add_mpi_exe_and_test(burgers1d_euler_tpetra ${ROOTNAME} burgers1d_euler_tpetra.cc 4 "PASSED")
ENDIF()
