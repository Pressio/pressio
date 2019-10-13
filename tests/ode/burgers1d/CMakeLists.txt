
add_subdirectory(eigen_explicit_euler)
add_subdirectory(eigen_implicit_euler)
add_subdirectory(eigen_implicit_euler_observer)
add_subdirectory(eigen_implicit_bdf2)
add_subdirectory(eigen_implicit_bdf2_observer)

# IF(HAVE_KOKKOS)
#   ADD_SUBDIRECTORY(kokkos_explicit_euler)
# ENDIF()

if(HAVE_TRILINOS)
  add_subdirectory(epetra_explicit_euler)
  add_subdirectory(tpetra_explicit_euler)
ENDIF()
