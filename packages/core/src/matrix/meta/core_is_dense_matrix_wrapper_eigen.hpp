
#ifndef CORE_IS_DENSE_MATRIX_WRAPPER_EIGEN_HPP_
#define CORE_IS_DENSE_MATRIX_WRAPPER_EIGEN_HPP_

#include "../core_matrix_traits.hpp"

namespace rompp{ namespace core{ namespace meta {

template <typename T, typename enable = void>
struct is_dense_matrix_wrapper_eigen : std::false_type {};

template <typename T>
struct is_dense_matrix_wrapper_eigen<
  T, ::rompp::mpl::enable_if_t<
       core::details::traits<T>::is_matrix &&
       (core::details::traits<T>::wrapped_matrix_identifier==
	core::details::WrappedMatrixIdentifier::DenseEigen)
       >
  >
  : std::true_type{};

}}}//end namespace rompp::core::meta
#endif
