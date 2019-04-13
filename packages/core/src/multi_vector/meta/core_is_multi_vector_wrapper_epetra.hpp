
#ifdef HAVE_TRILINOS
#ifndef CORE_IS_MULTI_VECTOR_WRAPPER_EPETRA_HPP_
#define CORE_IS_MULTI_VECTOR_WRAPPER_EPETRA_HPP_

#include "../core_multi_vector_traits.hpp"

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_multi_vector_wrapper_epetra : std::false_type {};

template <typename T>
struct is_multi_vector_wrapper_epetra<
  T, ::rompp::mpl::enable_if_t<
       core::details::traits<T>::is_multi_vector &&
       core::details::traits<T>::wrapped_multi_vector_identifier==
       core::details::WrappedMultiVectorIdentifier::Epetra
       >
  >
  : std::true_type{};

}}}//end namespace rompp::core::meta
#endif
#endif
