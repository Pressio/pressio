
#ifndef pressio_containers_multi_vector_include_HPP_
#define pressio_containers_multi_vector_include_HPP_

#include "./predicates/containers_native_eigen_multi_vector_meta.hpp"
#include "./predicates/containers_native_arbitrary_multi_vector_meta.hpp"
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "./predicates/containers_native_kokkos_multi_vector_meta.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./predicates/containers_native_epetra_multi_vector_meta.hpp"
#include "./predicates/containers_native_tpetra_block_multi_vector_meta.hpp"
#include "./predicates/containers_native_tpetra_multi_vector_meta.hpp"
#endif

#include "./predicates/containers_is_multi_vector_wrapper_eigen.hpp"
#include "./predicates/containers_is_multi_vector_wrapper_arbitrary.hpp"
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./predicates/containers_is_multi_vector_wrapper_epetra.hpp"
#include "./predicates/containers_is_multi_vector_wrapper_tpetra.hpp"
#include "./predicates/containers_is_multi_vector_wrapper_tpetra_block.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "./predicates/containers_is_multi_vector_wrapper_kokkos.hpp"
#endif
#include "./predicates/containers_is_multi_vector_wrapper.hpp"


#include "./containers_multi_vector_traits.hpp"


#include "./concrete/containers_multi_vector_arbitrary.hpp"
#include "./concrete/containers_multi_vector_sharedmem_eigen_dynamic.hpp"
#ifdef PRESSIO_ENABLE_TPL_TRILINOS
#include "./concrete/containers_multi_vector_distributed_epetra.hpp"
#include "./concrete/containers_multi_vector_distributed_tpetra_block.hpp"
#include "./concrete/containers_multi_vector_distributed_tpetra.hpp"
#endif
#ifdef PRESSIO_ENABLE_TPL_KOKKOS
#include "./concrete/containers_multi_vector_sharedmem_kokkos.hpp"
#endif

#endif


