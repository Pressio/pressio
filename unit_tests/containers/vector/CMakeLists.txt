
set(exeRootName containers_vector)

#############
### EIGEN ###
#############
if(PRESSIO_ENABLE_TPL_EIGEN)
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/eigen_dynamic_class.cc)
  add_serial_utest(${exeRootName}_eigen_dynamic_class ${SRC1})
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/eigen_static_class.cc)
  add_serial_utest(${exeRootName}_eigen_static_class ${SRC1})
endif()

##############
### KOKKOS ###
##############
if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_serial_utest_kokkos(${exeRootName}_kokkos_class
    ${CMAKE_CURRENT_SOURCE_DIR}/kokkos_class.cc )
endif()

################
### TRILINOS ###
################
if(PRESSIO_ENABLE_TPL_TRILINOS)
  ### TPETRA ###
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/tpetra_traits.cc)
  add_utest_mpi(${exeRootName}_tpetra_one ${SRC1} gTestMain_tpetra 3 )

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/tpetra_class_one.cc)
  add_utest_mpi(${exeRootName}_tpetra_two ${SRC1} gTestMain_tpetra 3 )

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/tpetra_class_two.cc)
  add_utest_mpi(${exeRootName}_tpetra_three ${SRC1} gTestMain_tpetra 3 )

  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/tpetra_manipulation.cc)
  add_utest_mpi(${exeRootName}_tpetra_manipulation ${SRC} gTestMain_tpetra 4)

  ### TPETRA BLOCK ###
  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/tpetra_block_traits.cc)
  add_utest_mpi(${exeRootName}_tpetra_block_one ${SRC1} gTestMain_tpetra 3 )

  set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/tpetra_block_class_two.cc)
  add_utest_mpi(${exeRootName}_tpetra_block_two ${SRC1} gTestMain_tpetra 3 )

  ### EPETRA ###
  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/epetra_traits.cc)
  add_serial_utest(${exeRootName}_epetra_traits ${SRC})

  set(SRC2 ${CMAKE_CURRENT_SOURCE_DIR}/epetra_class.cc)
  add_utest_mpi(${exeRootName}_epetra_class ${SRC2} gTestMain_mpi 3 )
endif()
