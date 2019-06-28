
#ifndef ALGEBRA_IS_MATRIX_WRAPPER_EIGEN_HPP_
#define ALGEBRA_IS_MATRIX_WRAPPER_EIGEN_HPP_

#include "algebra_is_dense_matrix_wrapper_eigen.hpp"
#include "algebra_is_sparse_matrix_wrapper_eigen.hpp"

namespace rompp{ namespace algebra{ namespace meta {

template <typename T, typename enable = void>
struct is_matrix_wrapper_eigen : std::false_type {};

template <typename T>
struct is_matrix_wrapper_eigen<
  T, ::rompp::mpl::enable_if_t<
       is_sparse_matrix_wrapper_eigen<T>::value or
       is_dense_matrix_wrapper_eigen<T>::value
       >
  >
  : std::true_type{};
//------------------------------------------------------------



}}}//end namespace rompp::algebra::meta
#endif
