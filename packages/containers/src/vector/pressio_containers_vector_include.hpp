
#ifndef pressio_containers_vector_include_HPP
#define pressio_containers_vector_include_HPP

#include "./predicates/containers_native_eigen_vector.hpp"
#include "./predicates/containers_native_arbitrary_vector.hpp"
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "./predicates/containers_native_kokkos_vector.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./predicates/containers_native_epetra_vector.hpp"
#include "./predicates/containers_native_teuchos_vector.hpp"
#include "./predicates/containers_native_tpetra_block_vector.hpp"
#include "./predicates/containers_native_tpetra_vector.hpp"
#endif

#include "./predicates/containers_is_vector_wrapper_arbitrary.hpp"
#include "./predicates/containers_is_vector_wrapper_eigen.hpp"
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./predicates/containers_is_dense_vector_wrapper_teuchos.hpp"
#include "./predicates/containers_is_vector_wrapper_epetra.hpp"
#include "./predicates/containers_is_vector_wrapper_tpetra_block.hpp"
#include "./predicates/containers_is_vector_wrapper_tpetra.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
#include "./predicates/containers_is_vector_wrapper_pybind.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "./predicates/containers_is_vector_wrapper_kokkos.hpp"
#endif
#include "./predicates/containers_is_vector_wrapper.hpp"

#include "./containers_vector_traits.hpp"

#include "./concrete/containers_vector_arbitrary.hpp"
#include "./concrete/containers_vector_sharedmem_eigen_dynamic.hpp"
#include "./concrete/containers_vector_sharedmem_eigen_static.hpp"
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
#include "./concrete/containers_vector_sharedmem_pybind11.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./concrete/containers_vector_distributed_epetra.hpp"
#include "./concrete/containers_vector_distributed_tpetra_block.hpp"
#include "./concrete/containers_vector_distributed_tpetra.hpp"
#include "./concrete/containers_vector_sharedmem_teuchos_serial_dense.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "./concrete/containers_vector_sharedmem_kokkos.hpp"
#endif

#endif
