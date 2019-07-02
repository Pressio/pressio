
#ifdef HAVE_ARMADILLO
#ifndef CONTAINERS_IS_VECTOR_WRAPPER_ARMADILLO_HPP_
#define CONTAINERS_IS_VECTOR_WRAPPER_ARMADILLO_HPP_

#include "../containers_vector_traits.hpp"

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_row_vector_wrapper_armadillo : std::false_type {};

template <typename T>
struct is_row_vector_wrapper_armadillo<
  T, ::pressio::mpl::enable_if_t<
       containers::details::traits<T>::is_vector &&
       containers::details::traits<T>::wrapped_vector_identifier==
       containers::details::WrappedVectorIdentifier::ArmadilloRow
       >
  > : std::true_type{};


template <typename T, typename enable = void>
struct is_column_vector_wrapper_armadillo : std::false_type {};

template <typename T>
struct is_column_vector_wrapper_armadillo<
  T, ::pressio::mpl::enable_if_t<
       containers::details::traits<T>::is_vector &&
       containers::details::traits<T>::wrapped_vector_identifier==
       containers::details::WrappedVectorIdentifier::ArmadilloCol
       >
  > : std::true_type{};

}}}//end namespace pressio::containers::meta
#endif
#endif // HAVE_ARMADILLO
