
#ifndef CORE_VECTOR_DOT_PRODUCT_HPP_
#define CORE_VECTOR_DOT_PRODUCT_HPP_

#include "core_vector_meta.hpp"
#include "core_vector_serial_eigen.hpp"
#include "core_vector_serial_stdlib.hpp"

namespace core{

template <typename vec_a_type,
	  typename vec_b_type,
	  typename res_t,
	  typename
	  std::enable_if<
	    // vec a has to be eigen 
	    details::traits<vec_a_type>::isVector &&
	    details::traits<vec_a_type>::isEigen && 
	    // vec b has to be eigen 
	    details::traits<vec_b_type>::isVector &&
	    details::traits<vec_b_type>::isEigen &&
	    // the two vectors need matching scalar type
	    std::is_same<typename
			 details::traits<vec_a_type>::scalar_t,
			 typename
			 details::traits<vec_b_type>::scalar_t
			 >::value &&
	    // result type must match too
	    std::is_same<typename
			 details::traits<vec_a_type>::scalar_t,
			 res_t>::value
	    >::type * = nullptr
	  >
void dotProduct(const vec_a_type & vecA,
		const vec_b_type & vecB,
		res_t & result)
{
  assert(vecA.size() == vecB.size());
  result = vecA.data()->dot(*vecB.data());
}

    
} // end namespace core
#endif

