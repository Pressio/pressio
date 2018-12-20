
#ifdef HAVE_TRILINOS
#ifndef CORE_NATIVE_TRILINOS_MATRIX_META_HPP_
#define CORE_NATIVE_TRILINOS_MATRIX_META_HPP_

#include "core_meta_basic.hpp"
#include "core_native_vector_meta.hpp"
#include <Epetra_CrsMatrix.h>
#include <Epetra_MultiVector.h>
#include <Epetra_LocalMap.h>
#include <Tpetra_CrsMatrix_decl.hpp>

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_matrix_sparse_distributed_tpetra : std::false_type {};

template <typename T>
struct is_matrix_sparse_distributed_tpetra<T,
      typename
      std::enable_if<
  std::is_same<T,
         Tpetra::CrsMatrix<
           typename T::impl_scalar_type,
           typename T::local_ordinal_type,
           typename T::global_ordinal_type,
           typename T::node_type
           >
         >::value
  >::type
      > : std::true_type{};
//--------------------------------------------

template <typename T, typename enable = void>
struct is_matrix_sparse_distributed_epetra
  : std::false_type {};

template<typename T>
struct is_matrix_sparse_distributed_epetra<T,
    typename std::enable_if<
      std::is_same<T, Epetra_CrsMatrix>::value
      >::type >
  : std::true_type{};

//----------------------------------------------------------------------

template <typename T, typename enable = void>
struct is_matrix_dense_distributed_epetra
  : std::false_type {};

template<typename T>
struct is_matrix_dense_distributed_epetra<T,
    typename std::enable_if<
      std::is_same<T, Epetra_MultiVector>::value
      >::type >
  : std::true_type{};


}}}//end namespace rompp::core::meta
#endif
#endif
