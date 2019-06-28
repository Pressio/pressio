
#ifndef ALGEBRA_IS_SPARSE_MATRIX_WRAPPER_EIGEN_HPP_
#define ALGEBRA_IS_SPARSE_MATRIX_WRAPPER_EIGEN_HPP_

#include "../algebra_matrix_traits.hpp"

namespace rompp{ namespace algebra{ namespace meta {

template <typename T, typename enable = void>
struct is_sparse_matrix_wrapper_eigen : std::false_type {};

template <typename T>
struct is_sparse_matrix_wrapper_eigen<
  T, ::rompp::mpl::enable_if_t<
       algebra::details::traits<T>::is_matrix &&
       (algebra::details::traits<T>::wrapped_matrix_identifier==
       algebra::details::WrappedMatrixIdentifier::SparseEigen)
       >
  >
  : std::true_type{};

}}}//end namespace rompp::algebra::meta
#endif
