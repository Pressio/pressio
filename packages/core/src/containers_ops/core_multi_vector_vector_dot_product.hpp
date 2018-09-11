
#ifndef CORE_MULTI_VECTOR_VECTOR_DOT_PRODUCT_HPP_
#define CORE_MULTI_VECTOR_VECTOR_DOT_PRODUCT_HPP_

#include "core_ops_meta.hpp"
#include "../vector/core_vector_meta.hpp"
#include "../multi_vector/core_multi_vector_meta.hpp"

namespace core{
namespace ops{

//  Epetra multivector with epetra vector
template <typename mvec_type,
	  typename vec_type,
    typename std::enable_if<
     core::meta::is_epetra_multi_vector_wrapper<mvec_type>::value &&
     core::meta::is_epetra_vector_wrapper<vec_type>::value &&
     core::meta::wrappers_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
void dot(const mvec_type & mvA,
	 const vec_type & vecB,
	 std::vector<typename
	 details::traits<mvec_type>::scalar_t> & result)
{
  ///computes dot product of each vector in mvA
  ///with vecB storing each value in result

  /* Apparently, trilinos does not support this... 
     the native Dot() method of multivectors is only for 
     dot product of two multivectors of the same size. 
     So we have to extract each column vector 
     from mvA and do dot product one a time*/
  
  // how many vectors are in mvA
  auto numVecs = mvA.globalNumVectors();
  if ( result.size() != numVecs )
    result.resize(numVecs);

  auto * mvNatData = mvA.data();
  const auto * vecNatData = vecB.data();
  for (decltype(numVecs) i=0; i<numVecs; i++){
    (*mvNatData)(i)->Dot(*vecNatData, &result[i]);
  }

}
    
//  Epetra multivector with epetra vector
template <typename mvec_type,
	  typename vec_type,
  typename std::enable_if<
    core::meta::is_epetra_multi_vector_wrapper<mvec_type>::value &&
    core::meta::is_epetra_vector_wrapper<vec_type>::value &&
    core::meta::wrappers_have_same_scalar<mvec_type, vec_type>::value
    > * = nullptr
  >
auto dot(const mvec_type & mvA,
	 const vec_type & vecB)
{
  using sc_t = typename details::traits<mvec_type>::scalar_t;
  // how many vectors are in mvA
  auto numVecs = mvA.globalNumVectors();
  using res_t = std::vector<sc_t>;
  res_t res(numVecs);
  dot(mvA, vecB, res);
  return res;
}

  
} // end namespace linalg
} // end namespace core
#endif
