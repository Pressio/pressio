
#ifndef CORE_EIGEN_DENSE_MAT_PRODUCT_EIGEN_DENSE_MAT_HPP_
#define CORE_EIGEN_DENSE_MAT_PRODUCT_EIGEN_DENSE_MAT_HPP_

#include "../core_ops_meta.hpp"
#include "../../matrix/core_matrix_meta.hpp"

namespace rompp{ namespace core{ namespace ops{

/*---------------------------------------------------------
C = A B
- A is dense matrix from eigen (can be static or dynamic)
- B is dense matrix from eigen (can be static or dynamic)
---------------------------------------------------------*/

template <typename T1, typename T2, typename T3,
core::meta::enable_if_t<
  core::meta::is_eigen_dense_matrix_wrapper<T1>::value &&
  core::meta::is_eigen_dense_matrix_wrapper<T2>::value &&
  core::meta::is_eigen_dense_matrix_wrapper<T3>::value &&
  core::meta::wrapper_triplet_have_same_scalar<T1,T2,T3>::value
  > * = nullptr
>
void product(const T1 & A, const T2 & B, T3 & C){

  assert(C.rows() == A.rows());
  assert(C.cols() == B.cols());
  (*C.data()) = (*A.data()) * (*B.data());
}


template <typename T1, typename T2,
     core::meta::enable_if_t<
       core::meta::is_eigen_dense_matrix_wrapper<T1>::value and
       core::meta::is_eigen_dense_matrix_wrapper<T2>::value and
       core::meta::wrapper_pair_have_same_scalar<T1,T2>::value
     > * = nullptr
    >
auto product(const T1 & A, const T2 & B)
-> core::Matrix<
    Eigen::Matrix<typename core::details::traits<T1>::scalar_t, -1, -1>
    >{

  using sc_t = typename core::details::traits<T1>::scalar_t;
  using nat_t = Eigen::Matrix<sc_t,-1,-1>;
  core::Matrix<nat_t> C(A.rows(),B.cols());
  product(A, B, C);
  return C;
}

}}}//end namespace rompp::core::ops
#endif
