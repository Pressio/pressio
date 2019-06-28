
#ifdef HAVE_BLAZE
#ifndef ALGEBRA_IS_VECTOR_WRAPPER_BLAZE_HPP_
#define ALGEBRA_IS_VECTOR_WRAPPER_BLAZE_HPP_

#include "../algebra_vector_traits.hpp"

namespace rompp{ namespace algebra{ namespace meta {

template <typename T, typename enable = void>
struct is_dynamic_vector_wrapper_blaze : std::false_type {};

template <typename T>
struct is_dynamic_vector_wrapper_blaze<
  T, ::rompp::mpl::enable_if_t<
       algebra::details::traits<T>::is_vector &&
       algebra::details::traits<T>::wrapped_vector_identifier==
       algebra::details::WrappedVectorIdentifier::BlazeDynamic
       >
  > : std::true_type{};


template <typename T, typename enable = void>
struct is_static_vector_wrapper_blaze : std::false_type {};

template <typename T>
struct is_static_vector_wrapper_blaze<
  T, ::rompp::mpl::enable_if_t<
       algebra::details::traits<T>::is_vector &&
       algebra::details::traits<T>::wrapped_vector_identifier==
       algebra::details::WrappedVectorIdentifier::BlazeStatic
       >
  > : std::true_type{};


}}}//end namespace rompp::algebra::meta
#endif
#endif
