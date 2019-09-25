#ifndef CONTAINERS_SRC_OPS_EIGEN_MAT_PROD_VECTOR_HPP_
#define CONTAINERS_SRC_OPS_EIGEN_MAT_PROD_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "containers_eigen_ops_helper_impl.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*
 * mat product vector
 */

/*---------------------------------------------
 * c = A * b
 * A : eigen sparse matrix wrapper
 * b: eigen vector wrapper
 * c: eigen vector
 *-----------------------------------------------*/
template <
  typename A_t,
  typename b_t,
  typename c_t,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_sparse_matrix_wrapper_eigen<A_t>::value &&
    containers::meta::is_vector_wrapper_eigen<b_t>::value &&
    containers::meta::is_vector_wrapper_eigen<c_t>::value &&
    containers::meta::wrapper_triplet_have_same_scalar<A_t,b_t,c_t>::value
    > * = nullptr
 >
void product(const A_t & A, const b_t & b, c_t & c){

  assert(A.cols() == b.size());
  assert(c.size() == A.rows());
  (*c.data()) = (*A.data()) * (*b.data());
}

/*---------------------------------------------
 * c = A * b
 * A : eigen sparse matrix wrapper
 * b: eigen vector wrapper
 * return c: eigen vector dynamic
 *-----------------------------------------------*/
template <
  typename A_t,
  typename b_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_sparse_matrix_wrapper_eigen<A_t>::value &&
    ::pressio::containers::meta::is_vector_wrapper_eigen<b_t>::value &&
    ::pressio::containers::details::traits<b_t>::is_dynamic and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<A_t, b_t>::value
    > * = nullptr
  >
b_t product(const A_t & A, const b_t & b)
{
  b_t c(A.rows());
  product(A,b,c);
  return c;
}



/*---------------------------------------------
 * c = A b
 * A : eigen dense matrix wrapper
 * b: eigen vector wrapper
 * c: eigen vector
 *-----------------------------------------------*/
template <
  typename A_t,
  typename b_t,
  typename c_t,
  bool transposeA = false,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_dense_matrix_wrapper_eigen<A_t>::value and
    containers::meta::is_vector_wrapper_eigen<b_t>::value and
    containers::meta::is_vector_wrapper_eigen<c_t>::value and
    containers::meta::wrapper_triplet_have_same_scalar<A_t, b_t, c_t>::value and
    transposeA == false
    > * = nullptr
  >
void product(const A_t & A, const b_t & b, c_t & c){

  assert(A.cols() == b.size());
  assert(c.size() == A.rows());
  (*c.data()) = (*A.data()) * (*b.data());
}


/*---------------------------------------------
 * c = A^T b
 * A : eigen dense matrix wrapper
 * b: eigen vector wrapper
 * c: eigen vector
 *-----------------------------------------------*/
template <
  typename A_t,
  typename b_t,
  typename c_t,
  bool transposeA = false,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_dense_matrix_wrapper_eigen<A_t>::value and
    containers::meta::is_vector_wrapper_eigen<b_t>::value and
    containers::meta::is_vector_wrapper_eigen<c_t>::value and
    containers::meta::wrapper_triplet_have_same_scalar<A_t, b_t, c_t>::value and
    transposeA == true
    > * = nullptr
  >
void product(const A_t & A, const b_t & b, c_t & c){

  assert(A.rows() == b.size());
  assert(c.size() == A.cols());
  (*c.data()) = (*A.data()).transpose() * (*b.data());
}


/*---------------------------------------------
 * c = A b
 * A : eigen dense matrix wrapper
 * b: eigen vector wrapper
 * c: eigen vector
 *-----------------------------------------------*/
template <
  typename A_t,
  typename b_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<A_t>::value and
    ::pressio::containers::meta::is_vector_wrapper_eigen<b_t>::value and
    ::pressio::containers::details::traits<b_t>::is_dynamic and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<A_t,b_t>::value
    > * = nullptr
   >
b_t product(const A_t & A, const b_t & b)
{
  b_t c(A.rows());
  product(A,b,c);
  return c;
}

}}}//end namespace pressio::containers::ops
#endif
