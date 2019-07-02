
#ifdef HAVE_TRILINOS
#ifndef CONTAINERS_TPETRA_BLOCK_MULTI_VECTOR_TPETRA_BLOCK_VECTOR_HPP_
#define CONTAINERS_TPETRA_BLOCK_MULTI_VECTOR_TPETRA_BLOCK_VECTOR_HPP_

#include "../../containers_ops_meta.hpp"
#include "../../../vector/containers_vector_meta.hpp"
#include "../../../multi_vector/containers_multi_vector_meta.hpp"

namespace pressio{ namespace containers{ namespace ops{

/*
 * overloads for c = A dot b
 * where:
 * A = wrapper Multivector
 * b = wrapper vector
 * c = a shared-mem vector wrapper
 */

//------------------------------------
// c = scalar *, passed in
//------------------------------------
template <typename mvec_type,
	  typename vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra_block<mvec_type>::value &&
    containers::meta::is_vector_wrapper_tpetra_block<vec_type>::value &&
    containers::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
void dot(const mvec_type & mvA,
	 const vec_type & vecB,
	 typename details::traits<mvec_type>::scalar_t * result){

  /* workaround the non-constness of getVectorView,
   * which is supposed to be const but it is not */
  using mv_tmp_t = Tpetra::Experimental::BlockVector<>;
  const auto vecB_vv = const_cast<mv_tmp_t*>(vecB.data())->getVectorView();

  const auto mvA_mvv = mvA.data()->getMultiVectorView();

  // how many vectors are in mvA
  const auto numVecs = mvA.globalNumVectors();
  for (auto i=0; i<numVecs; i++){
    // colI is a Teuchos::RCP<Vector<...>>
    const auto colI = mvA_mvv.getVector(i);
    result[i] = colI->dot(vecB_vv);
  }
}


//--------------------------------------------
// c = teuchos serial dense vector, passed in
//--------------------------------------------
template <typename mvec_type,
	  typename vec_type,
	  typename result_vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra_block<mvec_type>::value and
    containers::meta::is_vector_wrapper_tpetra_block<vec_type>::value and
    containers::meta::is_dense_vector_wrapper_teuchos<result_vec_type>::value and
    containers::meta::wrapper_triplet_have_same_scalar<mvec_type,
						 vec_type,
						 result_vec_type>::value and
    containers::details::traits<result_vec_type>::is_dynamic
    > * = nullptr
  >
void dot(const mvec_type & mvA,
	 const vec_type & vecB,
	 result_vec_type & result){
  const auto numVecs = mvA.globalNumVectors();

  if ( result.size() != numVecs )
    result.resize(numVecs);

  dot(mvA, vecB, result.data()->values());
}


//---------------------------------------
// c = eigen DYNAMIC vector, passed in
//---------------------------------------
template <typename mvec_type,
	  typename vec_type,
	  typename result_vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra_block<mvec_type>::value and
    containers::meta::is_vector_wrapper_tpetra_block<vec_type>::value and
    containers::meta::is_vector_wrapper_eigen<result_vec_type>::value and
    containers::meta::wrapper_triplet_have_same_scalar<mvec_type,
						 vec_type,
						 result_vec_type>::value and
    containers::details::traits<result_vec_type>::is_dynamic
    > * = nullptr
  >
void dot(const mvec_type & mvA,
	 const vec_type & vecB,
	 result_vec_type & result){

  ///computes dot product of each vector in mvA
  ///with vecB storing each value in result

  /* Apparently, trilinos does not support this...
     the native dot() method of multivectors is only for
     dot product of two multivectors with same # of columns.
     So we have to extract each column vector
     from mvA and do dot product one a time*/

  // how many vectors are in mvA
  auto numVecs = mvA.globalNumVectors();

  // check the result has right size
  if ( result.size() != numVecs )
    result.resize(numVecs);

  result.setZero();
  dot(mvA, vecB, result.data()->data());
}

//---------------------------------------
// c = eigen STATIC vector, passed in
//---------------------------------------
template <typename mvec_type,
	  typename vec_type,
	  typename result_vec_type,
  ::pressio::mpl::enable_if_t<
    containers::meta::is_multi_vector_wrapper_tpetra_block<mvec_type>::value and
    containers::meta::is_vector_wrapper_tpetra_block<vec_type>::value and
    containers::meta::is_vector_wrapper_eigen<result_vec_type>::value and
    containers::meta::wrapper_triplet_have_same_scalar<mvec_type,
						 vec_type,
						 result_vec_type>::value and
    containers::details::traits<result_vec_type>::is_static
    > * = nullptr
  >
void dot(const mvec_type & mvA,
	 const vec_type & vecB,
	 result_vec_type & result){

  ///computes dot product of each vector in mvA
  ///with vecB storing each value in result
  // check the result has right size
  assert( result.size() == mvA.globalNumVectors() );
  dot(mvA, vecB, result.data()->data());
}

}}} // end namespace pressio::containers::ops
#endif
#endif
