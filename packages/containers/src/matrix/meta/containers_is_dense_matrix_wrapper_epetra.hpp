
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_IS_DENSE_MATRIX_WRAPPER_EPETRA_HPP_
#define CONTAINERS_IS_DENSE_MATRIX_WRAPPER_EPETRA_HPP_

#include "../containers_matrix_traits.hpp"

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_dense_matrix_wrapper_epetra : std::false_type {};

template <typename T>
struct is_dense_matrix_wrapper_epetra<
  T, ::pressio::mpl::enable_if_t<
       containers::details::traits<T>::is_matrix &&
       containers::details::traits<T>::wrapped_matrix_identifier==
       containers::details::WrappedMatrixIdentifier::DenseEpetra
       >
  >
  : std::true_type{};

}}}//end namespace pressio::containers::meta
#endif
#endif
