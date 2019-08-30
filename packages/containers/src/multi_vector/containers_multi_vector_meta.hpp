
#ifndef CONTAINERS_MULTI_VECTOR_META_HPP_
#define CONTAINERS_MULTI_VECTOR_META_HPP_

#include "./meta/containers_is_multi_vector_wrapper.hpp"
#include "./meta/containers_is_multi_vector_wrapper_eigen.hpp"

#ifdef HAVE_TRILINOS
#include "./meta/containers_is_multi_vector_wrapper_epetra.hpp"
#include "./meta/containers_is_multi_vector_wrapper_tpetra.hpp"
#include "./meta/containers_is_multi_vector_wrapper_tpetra_block.hpp"
#endif

#ifdef HAVE_KOKKOS
#include "./meta/containers_is_multi_vector_wrapper_kokkos.hpp"
#endif

#endif
