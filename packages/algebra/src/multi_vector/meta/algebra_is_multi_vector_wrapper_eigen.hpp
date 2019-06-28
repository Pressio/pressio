
#ifndef ALGEBRA_IS_MULTI_VECTOR_WRAPPER_EIGEN_HPP_
#define ALGEBRA_IS_MULTI_VECTOR_WRAPPER_EIGEN_HPP_

#include "../algebra_multi_vector_traits.hpp"

namespace rompp{ namespace algebra{ namespace meta {

template <typename T, typename enable = void>
struct is_multi_vector_wrapper_eigen : std::false_type {};

template <typename T>
struct is_multi_vector_wrapper_eigen<
  T, ::rompp::mpl::enable_if_t<
       algebra::details::traits<T>::is_multi_vector &&
       algebra::details::traits<T>::wrapped_multi_vector_identifier==
       algebra::details::WrappedMultiVectorIdentifier::Eigen
       >
  >
  : std::true_type{};

}}}//end namespace rompp::algebra::meta
#endif
