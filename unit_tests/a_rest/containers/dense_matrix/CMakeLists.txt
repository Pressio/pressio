
set(exeRootName containers_dense_matrix)

#############
### EIGEN ###
#############
if(PRESSIO_ENABLE_TPL_EIGEN)
  # column-major matrix dynamic
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/eigen_dynamic.cc)
  add_serial_utest(${exeRootName}_eigen_dynamic_colmajor ${SRC1})
  target_compile_definitions(utest_${exeRootName}_eigen_dynamic_colmajor PRIVATE USE_COLMAJ)
  # row-major matrix dynamic
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/eigen_dynamic.cc)
  add_serial_utest(${exeRootName}_eigen_dynamic_rowmajor ${SRC1})
  target_compile_definitions(utest_${exeRootName}_eigen_dynamic_rowmajor PRIVATE USE_ROWMAJ)

  # column-major matrix static
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/eigen_static.cc)
  add_serial_utest(${exeRootName}_eigen_static_colmajor ${SRC1})
  target_compile_definitions(utest_${exeRootName}_eigen_static_colmajor PRIVATE USE_COLMAJ)
  # row-major matrix static
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/eigen_static.cc)
  add_serial_utest(${exeRootName}_eigen_static_rowmajor ${SRC1})
  target_compile_definitions(utest_${exeRootName}_eigen_static_rowmajor PRIVATE USE_ROWMAJ)
endif()

##############
### KOKKOS ###
##############
if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_serial_utest_kokkos(${exeRootName}_kokkos_strided_subview
    ${CMAKE_CURRENT_SOURCE_DIR}/kokkos_strided_subview.cc)

  add_serial_utest_kokkos(${exeRootName}_kokkos
    ${CMAKE_CURRENT_SOURCE_DIR}/kokkos_class.cc)
endif()

################
### TRILINOS ###
################
if(PRESSIO_ENABLE_TPL_TRILINOS)
  ### EPETRA ###
  set(SRC2 ${CMAKE_CURRENT_SOURCE_DIR}/epetra_class.cc)
  add_utest_mpi(${exeRootName}_epetra_class ${SRC2} gTestMain_mpi 3 )
endif()
