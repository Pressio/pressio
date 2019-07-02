
#ifdef HAVE_BLAZE
#ifndef CONTAINERS_IS_VECTOR_WRAPPER_BLAZE_HPP_
#define CONTAINERS_IS_VECTOR_WRAPPER_BLAZE_HPP_

#include "../containers_vector_traits.hpp"

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_dynamic_vector_wrapper_blaze : std::false_type {};

template <typename T>
struct is_dynamic_vector_wrapper_blaze<
  T, ::pressio::mpl::enable_if_t<
       containers::details::traits<T>::is_vector &&
       containers::details::traits<T>::wrapped_vector_identifier==
       containers::details::WrappedVectorIdentifier::BlazeDynamic
       >
  > : std::true_type{};


template <typename T, typename enable = void>
struct is_static_vector_wrapper_blaze : std::false_type {};

template <typename T>
struct is_static_vector_wrapper_blaze<
  T, ::pressio::mpl::enable_if_t<
       containers::details::traits<T>::is_vector &&
       containers::details::traits<T>::wrapped_vector_identifier==
       containers::details::WrappedVectorIdentifier::BlazeStatic
       >
  > : std::true_type{};


}}}//end namespace pressio::containers::meta
#endif
#endif
