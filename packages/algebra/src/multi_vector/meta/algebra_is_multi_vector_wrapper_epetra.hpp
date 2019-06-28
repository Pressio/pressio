
#ifdef HAVE_TRILINOS
#ifndef ALGEBRA_IS_MULTI_VECTOR_WRAPPER_EPETRA_HPP_
#define ALGEBRA_IS_MULTI_VECTOR_WRAPPER_EPETRA_HPP_

#include "../algebra_multi_vector_traits.hpp"

namespace rompp{ namespace algebra{ namespace meta {

template <typename T, typename enable = void>
struct is_multi_vector_wrapper_epetra : std::false_type {};

template <typename T>
struct is_multi_vector_wrapper_epetra<
  T, ::rompp::mpl::enable_if_t<
       algebra::details::traits<T>::is_multi_vector &&
       algebra::details::traits<T>::wrapped_multi_vector_identifier==
       algebra::details::WrappedMultiVectorIdentifier::Epetra
       >
  >
  : std::true_type{};

}}}//end namespace rompp::algebra::meta
#endif
#endif
