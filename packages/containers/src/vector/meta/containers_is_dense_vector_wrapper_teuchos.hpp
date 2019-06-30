
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_IS_DENSEWRAPPER_TEUCHOS_HPP_
#define CONTAINERS_IS_DENSEWRAPPER_TEUCHOS_HPP_

#include "../containers_vector_traits.hpp"

namespace rompp{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_dense_vector_wrapper_teuchos : std::false_type {};

template <typename T>
struct is_dense_vector_wrapper_teuchos<
  T, ::rompp::mpl::enable_if_t<
       containers::details::traits<T>::is_vector &&
       containers::details::traits<T>::wrapped_vector_identifier==
       containers::details::WrappedVectorIdentifier::TeuchosSerialDense
       >
  > : std::true_type{};

}}}//end namespace rompp::containers::meta
#endif
#endif
