#ifndef CONTAINERS_SRC_OPS_EIGEN_MAT_PROD_MULTI_VECTOR_HPP_
#define CONTAINERS_SRC_OPS_EIGEN_MAT_PROD_MULTI_VECTOR_HPP_

#include "../containers_ops_meta.hpp"
#include "containers_eigen_ops_helper_impl.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*
 * mat product multivector
 */

/*---------------------------------------------
 * C = A * B
 * A : eigen sparse matrix wrapper
 * B: eigen multivector wrapper
 * C: eigen dense matrix wrapper
 *-----------------------------------------------*/
template <
  typename mat_type,
  typename mvec_type,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_sparse_matrix_wrapper_eigen<mat_type>::value and
    ::pressio::containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value and
    ::pressio::containers::meta::wrapper_triplet_have_same_scalar<mat_type,
								  mvec_type,
								  result_t>::value
    > * = nullptr
  >
void product(const mat_type & A, const mvec_type & mv, result_t & C)
{
  assert( C.data()->rows() == A.data()->rows() );
  assert( mv.length() == A.cols() );
  assert( C.data()->cols() == mv.numVectors() );
  (*C.data()) = (*A.data()) * (*mv.data());
}//end function



/*---------------------------------------------
 * C = A * B
 * A : eigen sparse matrix wrapper
 * B: eigen multivector wrapper
 * C: eigen multivector wrapper
 *-----------------------------------------------*/
template <
  typename mat_type,
  typename mvec_type,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_sparse_matrix_wrapper_eigen<mat_type>::value and
    ::pressio::containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mat_type, mvec_type>::value
    > * = nullptr
  >
void product(const mat_type & A, const mvec_type & mv, mvec_type & C)
{
  assert( C.data()->rows() == A.data()->rows() );
  assert( mv.length() == A.cols() );
  assert( C.data()->cols() == mv.numVectors() );
  (*C.data()) = (*A.data()) * (*mv.data());
}//end function



/*---------------------------------------------
 * C = A * B
 * A : eigen sparse matrix wrapper
 * B: eigen multivector wrapper
 * C: is returned and is dynamic
 *-----------------------------------------------*/
template <
  typename mat_type,
  typename mvec_type,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_sparse_matrix_wrapper_eigen<mat_type>::value and
    ::pressio::containers::meta::is_multi_vector_wrapper_eigen<mvec_type>::value and
    (::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value or
     ::pressio::containers::meta::is_multi_vector_wrapper_eigen<result_t>::value) and
     ::pressio::containers::details::traits<result_t>::is_dynamic and
     ::pressio::containers::meta::wrapper_triplet_have_same_scalar<mat_type, mvec_type, result_t>::value
    > * = nullptr
  >
result_t product(const mat_type & A, const mvec_type & mv)
{
  result_t C(A.rows(), mv.numVectors());
  product(A,mv,C);
  return C;
}//end function


}}}//end namespace pressio::containers::ops
#endif
