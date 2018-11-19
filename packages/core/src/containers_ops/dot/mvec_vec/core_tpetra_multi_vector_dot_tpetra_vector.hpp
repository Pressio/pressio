
#ifndef CORE_TPETRA_MULTI_VECTOR_TPETRA_VECTOR_HPP_
#define CORE_TPETRA_MULTI_VECTOR_TPETRA_VECTOR_HPP_

#include "../../core_ops_meta.hpp"
#include "../../../vector/core_vector_meta.hpp"
#include "../../../multi_vector/core_multi_vector_meta.hpp"

namespace rompp{ namespace core{ namespace ops{

// Tpetra multivector dot tpetra vector
// result stored into Eigen vector passed by reference

#ifdef HAVE_TRILINOS
template <typename mvec_type,
	  typename vec_type,
	  typename result_vec_type,
  core::meta::enable_if_t<
    core::meta::is_tpetra_multi_vector_wrapper<mvec_type>::value and
    core::meta::is_tpetra_vector_wrapper<vec_type>::value and
    core::meta::wrapper_pair_have_same_scalar<mvec_type, vec_type>::value and
    core::meta::is_eigen_vector_wrapper<result_vec_type>::value and
    core::meta::wrapper_pair_have_same_scalar<vec_type,result_vec_type>::value
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

  for (decltype(numVecs) i=0; i<numVecs; i++){
    // colI is a Teuchos::RCP<Vector<...>>
    auto colI = mvA.data()->getVector(i);
    result[i] = colI->dot(*vecB.data());
  }
}
#endif
//--------------------------------------------------------


}}} // end namespace rompp::core::ops
#endif
