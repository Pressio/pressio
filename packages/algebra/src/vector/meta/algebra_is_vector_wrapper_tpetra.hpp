
#ifdef HAVE_TRILINOS
#ifndef ALGEBRA_IS_VECTOR_WRAPPER_TPETRA_HPP_
#define ALGEBRA_IS_VECTOR_WRAPPER_TPETRA_HPP_

#include "../algebra_vector_traits.hpp"

namespace rompp{ namespace algebra{ namespace meta {

template <typename T, typename enable = void>
struct is_vector_wrapper_tpetra : std::false_type {};

template <typename T>
struct is_vector_wrapper_tpetra<
  T, ::rompp::mpl::enable_if_t<
       algebra::details::traits<T>::is_vector &&
       algebra::details::traits<T>::wrapped_vector_identifier==
       algebra::details::WrappedVectorIdentifier::Tpetra
       >
  > : std::true_type{};

}}}//end namespace rompp::algebra::meta
#endif
#endif
