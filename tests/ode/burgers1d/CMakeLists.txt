
if(PRESSIO_ENABLE_TPL_EIGEN)
  add_subdirectory(eigen_explicit_euler)
  add_subdirectory(eigen_implicit_euler)
  add_subdirectory(eigen_implicit_euler_observer)
  add_subdirectory(eigen_implicit_bdf2)
  add_subdirectory(eigen_implicit_bdf2_observer)
  add_subdirectory(eigen_implicit_crank_nicolson)
endif()

if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_subdirectory(kokkos_explicit_euler)
endif()

if(PRESSIO_ENABLE_TPL_TRILINOS)
  add_subdirectory(epetra_explicit_euler)
  # add_subdirectory(tpetra_explicit_euler)
ENDIF()
