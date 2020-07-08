
#ifndef pressio_containers_matrix_include_HPP
#define pressio_containers_matrix_include_HPP

#include "./predicates/containers_native_eigen_matrix_meta.hpp"
#include "./predicates/containers_native_kokkos_matrix_meta.hpp"
#include "./predicates/containers_native_trilinos_matrix_meta.hpp"
#include "./predicates/containers_native_arbitrary_matrix_meta.hpp"

#include "./predicates/containers_is_dense_matrix_wrapper_eigen.hpp"
#include "./predicates/containers_is_dense_matrix_wrapper_epetra.hpp"
#include "./predicates/containers_is_dense_matrix_wrapper_kokkos.hpp"
#include "./predicates/containers_is_dense_matrix_wrapper_teuchos.hpp"
#include "./predicates/containers_is_sparse_matrix_wrapper_eigen.hpp"
#include "./predicates/containers_is_sparse_matrix_wrapper_epetra.hpp"
#include "./predicates/containers_is_sparse_matrix_wrapper_kokkos.hpp"
#include "./predicates/containers_is_sparse_matrix_wrapper_tpetra.hpp"
#include "./predicates/containers_is_matrix_wrapper_arbitrary.hpp"
#include "./predicates/containers_is_matrix_wrapper_eigen.hpp"
#include "./predicates/containers_is_matrix_wrapper_kokkos.hpp"
#include "./predicates/containers_is_matrix_wrapper_pybind.hpp"
#include "./predicates/containers_is_matrix_wrapper.hpp"

#include "./containers_matrix_traits.hpp"

#include "./concrete/containers_matrix_arbitrary.hpp"
#include "./concrete/containers_matrix_sharedmem_pybind11.hpp"
#include "./concrete/containers_matrix_dense_distributed_epetra.hpp"
#include "./concrete/containers_matrix_dense_sharedmem_eigen_dynamic.hpp"
#include "./concrete/containers_matrix_dense_sharedmem_eigen_static.hpp"
#include "./concrete/containers_matrix_dense_sharedmem_kokkos.hpp"
#include "./concrete/containers_matrix_dense_sharedmem_teuchos_serial.hpp"
#include "./concrete/containers_matrix_sparse_sharedmem_eigen.hpp"

#endif
