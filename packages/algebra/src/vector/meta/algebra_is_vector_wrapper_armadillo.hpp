
#ifdef HAVE_ARMADILLO
#ifndef ALGEBRA_IS_VECTOR_WRAPPER_ARMADILLO_HPP_
#define ALGEBRA_IS_VECTOR_WRAPPER_ARMADILLO_HPP_

#include "../algebra_vector_traits.hpp"

namespace rompp{ namespace algebra{ namespace meta {

template <typename T, typename enable = void>
struct is_row_vector_wrapper_armadillo : std::false_type {};

template <typename T>
struct is_row_vector_wrapper_armadillo<
  T, ::rompp::mpl::enable_if_t<
       algebra::details::traits<T>::is_vector &&
       algebra::details::traits<T>::wrapped_vector_identifier==
       algebra::details::WrappedVectorIdentifier::ArmadilloRow
       >
  > : std::true_type{};


template <typename T, typename enable = void>
struct is_column_vector_wrapper_armadillo : std::false_type {};

template <typename T>
struct is_column_vector_wrapper_armadillo<
  T, ::rompp::mpl::enable_if_t<
       algebra::details::traits<T>::is_vector &&
       algebra::details::traits<T>::wrapped_vector_identifier==
       algebra::details::WrappedVectorIdentifier::ArmadilloCol
       >
  > : std::true_type{};

}}}//end namespace rompp::algebra::meta
#endif
#endif // HAVE_ARMADILLO
