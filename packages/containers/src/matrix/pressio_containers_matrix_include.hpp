
#ifndef pressio_containers_matrix_include_HPP
#define pressio_containers_matrix_include_HPP

// predicates native
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "./predicates/containers_native_kokkos_matrix_meta.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./predicates/containers_native_trilinos_matrix_meta.hpp"
#endif
#include "./predicates/containers_native_eigen_matrix_meta.hpp"
#include "./predicates/containers_native_arbitrary_matrix_meta.hpp"

// predicates wrappers
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./predicates/containers_is_dense_matrix_wrapper_epetra.hpp"
#include "./predicates/containers_is_dense_matrix_wrapper_teuchos.hpp"
#include "./predicates/containers_is_sparse_matrix_wrapper_epetra.hpp"
#include "./predicates/containers_is_sparse_matrix_wrapper_tpetra.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "./predicates/containers_is_dense_matrix_wrapper_kokkos.hpp"
#include "./predicates/containers_is_sparse_matrix_wrapper_kokkos.hpp"
#include "./predicates/containers_is_matrix_wrapper_kokkos.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
#include "./predicates/containers_is_matrix_wrapper_pybind11.hpp"
#endif
#include "./predicates/containers_is_sparse_matrix_wrapper_eigen.hpp"
#include "./predicates/containers_is_dense_matrix_wrapper_eigen.hpp"
#include "./predicates/containers_is_matrix_wrapper_eigen.hpp"
#include "./predicates/containers_is_matrix_wrapper_arbitrary.hpp"
#include "./predicates/containers_is_matrix_wrapper.hpp"

// traits
#include "./containers_matrix_traits.hpp"

// concrete
#ifdef PRESSIO_ENABLE_TPL_PYBIND11
#include "./concrete/containers_matrix_sharedmem_pybind11.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./concrete/containers_matrix_dense_distributed_epetra.hpp"
#include "./concrete/containers_matrix_dense_sharedmem_teuchos_serial.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "./concrete/containers_matrix_dense_sharedmem_kokkos.hpp"
#endif
#include "./concrete/containers_matrix_arbitrary.hpp"
#include "./concrete/containers_matrix_sparse_sharedmem_eigen.hpp"
#include "./concrete/containers_matrix_dense_sharedmem_eigen_dynamic.hpp"
#include "./concrete/containers_matrix_dense_sharedmem_eigen_static.hpp"

#endif
