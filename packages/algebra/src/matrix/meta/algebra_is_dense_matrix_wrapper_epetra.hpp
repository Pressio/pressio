
#ifdef HAVE_TRILINOS
#ifndef ALGEBRA_IS_DENSE_MATRIX_WRAPPER_EPETRA_HPP_
#define ALGEBRA_IS_DENSE_MATRIX_WRAPPER_EPETRA_HPP_

#include "../algebra_matrix_traits.hpp"

namespace rompp{ namespace algebra{ namespace meta {

template <typename T, typename enable = void>
struct is_dense_matrix_wrapper_epetra : std::false_type {};

template <typename T>
struct is_dense_matrix_wrapper_epetra<
  T, ::rompp::mpl::enable_if_t<
       algebra::details::traits<T>::is_matrix &&
       algebra::details::traits<T>::wrapped_matrix_identifier==
       algebra::details::WrappedMatrixIdentifier::DenseEpetra
       >
  >
  : std::true_type{};

}}}//end namespace rompp::algebra::meta
#endif
#endif
