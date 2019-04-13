
#ifdef HAVE_ARMADILLO
#ifndef CORE_IS_VECTOR_WRAPPER_ARMADILLO_HPP_
#define CORE_IS_VECTOR_WRAPPER_ARMADILLO_HPP_

#include "../core_vector_traits.hpp"

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_row_vector_wrapper_armadillo : std::false_type {};

template <typename T>
struct is_row_vector_wrapper_armadillo<
  T, ::rompp::mpl::enable_if_t<
       core::details::traits<T>::is_vector &&
       core::details::traits<T>::wrapped_vector_identifier==
       core::details::WrappedVectorIdentifier::ArmadilloRow
       >
  > : std::true_type{};


template <typename T, typename enable = void>
struct is_column_vector_wrapper_armadillo : std::false_type {};

template <typename T>
struct is_column_vector_wrapper_armadillo<
  T, ::rompp::mpl::enable_if_t<
       core::details::traits<T>::is_vector &&
       core::details::traits<T>::wrapped_vector_identifier==
       core::details::WrappedVectorIdentifier::ArmadilloCol
       >
  > : std::true_type{};

}}}//end namespace rompp::core::meta
#endif
#endif // HAVE_ARMADILLO
