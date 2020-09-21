

set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_serial_eigen_dynamic_class.cc)
add_serial_utest(containers_vector_serial "${SRC1}")

set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_sharedmem_eigen_dynamic_constructors_and_smf.cc)
add_serial_utest(containers_vector_sharedmem_eigen_dynamic_constructors_and_smf ${SRC1})

set(SRC1 ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_sharedmem_eigen_static_constructors_and_smf.cc)
add_serial_utest(containers_vector_sharedmem_eigen_static_constructors_and_smf ${SRC1})



if(PRESSIO_ENABLE_TPL_PYBIND11)
  set(SRC2 ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_pybind.cc)
  add_serial_utest(containers_vector_pybind11 "${SRC2}")
endif()


if(PRESSIO_ENABLE_TPL_KOKKOS)
  add_serial_utest_kokkos(containers_vector_sharedmem_kokkos_constructors_and_smf
    containers_vector_sharedmem_kokkos_constructors_and_smf.cc )

  add_serial_utest_kokkos(containers_vector_sharedmem_kokkos_class
    containers_vector_sharedmem_kokkos_class.cc )

  add_serial_utest_kokkos(containers_vector_sharedmem_kokkos_static_dynamic
    containers_vector_sharedmem_kokkos_static_dynamic.cc)
endif()


if(PRESSIO_ENABLE_TPL_TRILINOS)
  set(SRC1
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_distributed_tpetra_traits.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_distributed_tpetra_class.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_distributed_tpetra_block_traits.cc
    ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_distributed_tpetra_block_class.cc
    )
  add_utest_mpi(containers_vector_tpetra "${SRC1}" gTestMain_tpetra 3 )

  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_manipulation_tpetra.cc)
  add_utest_mpi(containers_vector_manipulation_tpetra ${SRC} gTestMain_tpetra 4)

  set(SRC2 ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_distributed_tpetra_constructors_and_smf.cc)
  add_utest_mpi(containers_vector_tpetra_constr_and_smf "${SRC2}" gTestMain_tpetra 3 )

  set(SRC ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_distributed_epetra_traits.cc)
  add_serial_utest(containers_vector_distributed_epetra_traits ${SRC})

  set(SRC2 ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_distributed_epetra_class.cc)
  add_utest_mpi(containers_vector_epetra "${SRC2}" gTestMain_mpi 3 )

  set(SRC2 ${CMAKE_CURRENT_SOURCE_DIR}/containers_vector_distributed_epetra_constructors_and_smf.cc)
  add_utest_mpi(containers_vector_epetra_constr_and_smf "${SRC2}" gTestMain_mpi 3 )
endif()
