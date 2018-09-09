
#ifndef CORE_MULTI_VECTOR_VECTOR_DOT_PRODUCT_HPP_
#define CORE_MULTI_VECTOR_VECTOR_DOT_PRODUCT_HPP_

#include "../core_multi_vector_traits.hpp"
#include "../../meta/core_vector_meta.hpp"
#include "../../meta/core_multi_vector_meta.hpp"

namespace core{
namespace multivec_ops{

//  Epetra multivector with epetra vector
template <typename mvec_type,
	  typename vec_type,
	  core::meta::enable_if_t<
	   details::traits<mvec_type>::isMultiVector &&
	   details::traits<vec_type>::isEpetra &&
	   details::traits<vec_type>::isVector &&
	   details::traits<vec_type>::isEpetra &&
	   std::is_same<
	     typename details::traits<mvec_type>::scalar_t,
	     typename details::traits<vec_type>::scalar_t
	     >::value
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
	  core::meta::enable_if_t<
	   details::traits<mvec_type>::isMultiVector &&
	   details::traits<vec_type>::isEpetra && 
	   details::traits<vec_type>::isVector &&
	   details::traits<vec_type>::isEpetra &&
	   std::is_same<
	     typename details::traits<mvec_type>::scalar_t,
	     typename details::traits<vec_type>::scalar_t
	     >::value 
	   > * = nullptr
	  >
auto dot(const mvec_type & mvA,
	 const vec_type & vecB)
{
  using sc_t = typename details::traits<mvec_type>::scalar_t;
  // how many vectors are in mvA
  auto numVecs = mvA.globalNumVectors();
  std::vector<sc_t> res(numVecs);
  dot(mvA, vecB, res);
  return res;
}

  
} // end namespace multivec_ops
} // end namespace core
#endif
