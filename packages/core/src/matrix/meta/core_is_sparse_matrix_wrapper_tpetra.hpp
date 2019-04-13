
#ifdef HAVE_TRILINOS
#ifndef CORE_IS_SPARSE_MATRIX_WRAPPER_TPETRA_HPP_
#define CORE_IS_SPARSE_MATRIX_WRAPPER_TPETRA_HPP_

#include "../core_matrix_traits.hpp"

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_sparse_matrix_wrapper_tpetra : std::false_type {};

template <typename T>
struct is_sparse_matrix_wrapper_tpetra<
  T, ::rompp::mpl::enable_if_t<
       core::details::traits<T>::is_matrix &&
       core::details::traits<T>::wrapped_matrix_identifier==
       core::details::WrappedMatrixIdentifier::SparseTpetra
       >
  >
  : std::true_type{};

}}}//end namespace rompp::core::meta
#endif
#endif
