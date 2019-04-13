
#ifdef HAVE_TRILINOS
#ifndef CORE_IS_DENSE_MATRIX_WRAPPER_EPETRA_HPP_
#define CORE_IS_DENSE_MATRIX_WRAPPER_EPETRA_HPP_

#include "../core_matrix_traits.hpp"

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_dense_matrix_wrapper_epetra : std::false_type {};

template <typename T>
struct is_dense_matrix_wrapper_epetra<
  T, ::rompp::mpl::enable_if_t<
       core::details::traits<T>::is_matrix &&
       core::details::traits<T>::wrapped_matrix_identifier==
       core::details::WrappedMatrixIdentifier::DenseEpetra
       >
  >
  : std::true_type{};

}}}//end namespace rompp::core::meta
#endif
#endif
