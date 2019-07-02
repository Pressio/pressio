
#ifndef CONTAINERS_EIGEN_MAT_PRODUCT_MAT_HPP_
#define CONTAINERS_EIGEN_MAT_PRODUCT_MAT_HPP_

#include "../containers_ops_meta.hpp"
#include "../../matrix/containers_matrix_meta.hpp"
#include "../containers_eigen_ops_helper_impl.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*---------------------------------------------------------
 *
 * C = A B
 * A matrix wrapper for an eigen matrix (static/dynamic)
 * B matrix wrapper for an eigen matrix (static/dynamic)
 *
 *--------------------------------------------------------*/


template <
  typename TA, typename TB, typename TC,
  bool transposeA = false, bool transposeB = false,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_matrix_wrapper_eigen<TA>::value and
    containers::meta::is_matrix_wrapper_eigen<TB>::value and
    containers::meta::is_matrix_wrapper_eigen<TC>::value and
    containers::meta::wrapper_triplet_have_same_scalar<TA,TB,TC>::value
    > * = nullptr
  >
void product(const TA & A, const TB & B, TC & C){

  using implClass_t = impl::eig_mat_mat_product<transposeA,
						transposeB>;
  implClass_t()(A,B,C);
}

// DENSE times DENSE
template <
  typename T1, typename T2,
    bool transposeA = false, bool transposeB = false,
    ::pressio::mpl::enable_if_t<
      containers::meta::is_dense_matrix_wrapper_eigen<T1>::value and
      containers::meta::is_dense_matrix_wrapper_eigen<T2>::value and
      containers::meta::wrapper_pair_have_same_scalar<T1,T2>::value
      > * = nullptr
    >
auto product(const T1 & A, const T2 & B)
  -> typename impl::eigenMatMatProdRetTypeHelper<T1, T2>::prod_type
{
  using implClass_t = impl::eig_mat_mat_product<transposeA,
						transposeB>;
  return implClass_t()(A,B);
}

// SPARSE times DENSE
template <
  typename T1, typename T2,
    bool transposeA = false, bool transposeB = false,
    ::pressio::mpl::enable_if_t<
      containers::meta::is_sparse_matrix_wrapper_eigen<T1>::value and
      containers::meta::is_dense_matrix_wrapper_eigen<T2>::value and
      containers::meta::wrapper_pair_have_same_scalar<T1,T2>::value
      > * = nullptr
    >
auto product(const T1 & A, const T2 & B)
  -> typename impl::eigenMatMatProdRetTypeHelper<T1, T2>::prod_type
{
  using implClass_t = impl::eig_mat_mat_product<transposeA,
						transposeB>;
  return implClass_t()(A,B);
}

// DENSE times SPARSE
template <
  typename T1, typename T2,
    bool transposeA = false, bool transposeB = false,
    ::pressio::mpl::enable_if_t<
      containers::meta::is_dense_matrix_wrapper_eigen<T1>::value and
      containers::meta::is_sparse_matrix_wrapper_eigen<T2>::value and
      containers::meta::wrapper_pair_have_same_scalar<T1,T2>::value
      > * = nullptr
    >
auto product(const T1 & A, const T2 & B)
  -> typename impl::eigenMatMatProdRetTypeHelper<T1, T2>::prod_type
{
  using implClass_t = impl::eig_mat_mat_product<transposeA,
						transposeB>;
  return implClass_t()(A,B);
}


// SPARSE times SPARSE
template <
  typename T1, typename T2,
    bool transposeA = false, bool transposeB = false,
    ::pressio::mpl::enable_if_t<
      containers::meta::is_sparse_matrix_wrapper_eigen<T1>::value and
      containers::meta::is_sparse_matrix_wrapper_eigen<T2>::value and
      containers::meta::wrapper_pair_have_same_scalar<T1,T2>::value
      > * = nullptr
    >
auto product(const T1 & A, const T2 & B)
  -> typename impl::eigenMatMatProdRetTypeHelper<T1, T2>::prod_type
{
  using implClass_t = impl::eig_mat_mat_product<transposeA,
						transposeB>;
  return implClass_t()(A,B);
}


}}} // end namespace pressio::containers::ops
#endif
