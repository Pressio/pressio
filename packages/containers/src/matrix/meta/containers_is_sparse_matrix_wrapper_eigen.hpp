
#ifndef CONTAINERS_IS_SPARSE_MATRIX_WRAPPER_EIGEN_HPP_
#define CONTAINERS_IS_SPARSE_MATRIX_WRAPPER_EIGEN_HPP_

#include "../containers_matrix_traits.hpp"

namespace pressio{ namespace containers{ namespace meta {

template <typename T, typename enable = void>
struct is_sparse_matrix_wrapper_eigen : std::false_type {};

template <typename T>
struct is_sparse_matrix_wrapper_eigen<
  T, ::pressio::mpl::enable_if_t<
       containers::details::traits<T>::is_matrix &&
       (containers::details::traits<T>::wrapped_matrix_identifier==
       containers::details::WrappedMatrixIdentifier::SparseEigen)
       >
  >
  : std::true_type{};

}}}//end namespace pressio::containers::meta
#endif
