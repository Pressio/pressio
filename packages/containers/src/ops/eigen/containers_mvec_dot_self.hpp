#ifndef CONTAINERS_SRC_OPS_EIGEN_MULTI_VECTOR_DOT_SELF_HPP_
#define CONTAINERS_SRC_OPS_EIGEN_MULTI_VECTOR_DOT_SELF_HPP_

#include "../containers_ops_meta.hpp"
#include "../../multi_vector/containers_multi_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*
 * multi_vector dot self
 */

/* void specialize for:
 * result_t = dense dynamic eigen matrix wrapper
 */
template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_eigen<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value and
    ::pressio::containers::details::traits<result_t>::is_dynamic and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mvec_t, result_t>::value
    > * = nullptr
  >
void dot_self(const mvec_t & A, result_t & C)
{
  const auto nAcols = A.data()->cols();
  // since C is dynamic, I can resize if needed
  if(C.data()->rows() != nAcols || C.data()->cols() != nAcols)
    C.data()->resize( nAcols, nAcols );

  *C.data() = A.data()->transpose() * (*A.data());
}


/* non void specialize for:
 * result_t = dense dynamic eigen matrix wrapper
 */
template <
  typename mvec_t,
  typename result_t,
  ::pressio::mpl::enable_if_t<
    ::pressio::containers::meta::is_multi_vector_wrapper_eigen<mvec_t>::value and
    ::pressio::containers::meta::is_dense_matrix_wrapper_eigen<result_t>::value and
    ::pressio::containers::details::traits<result_t>::is_dynamic and
    ::pressio::containers::meta::wrapper_pair_have_same_scalar<mvec_t, result_t>::value
    > * = nullptr
  >
result_t dot_self(const mvec_t & A)
{
  const auto numVecsA = A.numVectors();
  result_t C(numVecsA, numVecsA);
  dot_self(A, C);
  return C;
}


}}}//end namespace pressio::containers::ops
#endif
