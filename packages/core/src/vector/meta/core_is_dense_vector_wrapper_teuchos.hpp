
#ifdef HAVE_TRILINOS
#ifndef CORE_IS_DENSEWRAPPER_TEUCHOS_HPP_
#define CORE_IS_DENSEWRAPPER_TEUCHOS_HPP_

#include "../core_vector_traits.hpp"

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_dense_vector_wrapper_teuchos : std::false_type {};

template <typename T>
struct is_dense_vector_wrapper_teuchos<
  T, ::rompp::mpl::enable_if_t<
       core::details::traits<T>::is_vector &&
       core::details::traits<T>::wrapped_vector_identifier==
       core::details::WrappedVectorIdentifier::TeuchosSerialDense
       >
  > : std::true_type{};

}}}//end namespace rompp::core::meta
#endif
#endif
