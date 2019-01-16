
#ifndef CORE_MATRIX_OPERATIONS_SPARSE_SPARSE_MATRIX_PRODUCT_HPP_
#define CORE_MATRIX_OPERATIONS_SPARSE_SPARSE_MATRIX_PRODUCT_HPP_

#include "../core_ops_meta.hpp"
#include "../../matrix/core_matrix_meta.hpp"
#include "core_eigen_mat_mat_impl.hpp"

namespace rompp{ namespace core{ namespace ops{

/*---------------------------------------------------------
C = A B
- A is sparse matrix from eigen
- B is sparse matrix from eigen
---------------------------------------------------------*/

template <
  typename TA, typename TB, typename TC,
  bool transposeA = false, bool transposeB = false,
  core::meta::enable_if_t<
    core::meta::is_eigen_sparse_matrix_wrapper<TA>::value &&
    core::meta::is_eigen_sparse_matrix_wrapper<TB>::value &&
    core::meta::is_eigen_sparse_matrix_wrapper<TC>::value &&
    core::meta::wrapper_triplet_have_same_scalar<TA,TB,TC>::value
    > * = nullptr
  >
void product(const TA & A, const TB & B, TC & C){

  using implClass_t = impl::eig_mat_mat_product<transposeA,
						transposeB>;
  implClass_t()(A,B,C);
  // assert(C.rows() == A.rows());
  // assert(C.cols() == B.cols());
  // (*C.data()) = (*A.data()) * (*B.data());
}


template <
  typename TA, typename TB,
    bool transposeA = false, bool transposeB = false,
    core::meta::enable_if_t<
      core::meta::is_eigen_sparse_matrix_wrapper<TA>::value &&
      core::meta::is_eigen_sparse_matrix_wrapper<TB>::value &&
      core::meta::wrapper_pair_have_same_scalar<TA,TB>::value
      > * = nullptr
    >
TA product(const TA & A, const TB & B){

  using implClass_t = impl::eig_mat_mat_product<transposeA,
						transposeB>;
  return implClass_t()(A,B);

  // TA C(A.rows(),B.cols());
  // product(A, B, C);
  // return C;
}


}}} // end namespace rompp::core::ops
#endif
