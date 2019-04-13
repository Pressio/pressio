
#ifdef HAVE_BLAZE
#ifndef CORE_IS_VECTOR_WRAPPER_BLAZE_HPP_
#define CORE_IS_VECTOR_WRAPPER_BLAZE_HPP_

#include "../core_vector_traits.hpp"

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_dynamic_vector_wrapper_blaze : std::false_type {};

template <typename T>
struct is_dynamic_vector_wrapper_blaze<
  T, ::rompp::mpl::enable_if_t<
       core::details::traits<T>::is_vector &&
       core::details::traits<T>::wrapped_vector_identifier==
       core::details::WrappedVectorIdentifier::BlazeDynamic
       >
  > : std::true_type{};


template <typename T, typename enable = void>
struct is_static_vector_wrapper_blaze : std::false_type {};

template <typename T>
struct is_static_vector_wrapper_blaze<
  T, ::rompp::mpl::enable_if_t<
       core::details::traits<T>::is_vector &&
       core::details::traits<T>::wrapped_vector_identifier==
       core::details::WrappedVectorIdentifier::BlazeStatic
       >
  > : std::true_type{};


}}}//end namespace rompp::core::meta
#endif
#endif
