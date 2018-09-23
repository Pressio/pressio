
#ifndef CORE_VECTOR_OPERATIONS_VECTOR_DOT_PRODUCT_HPP_
#define CORE_VECTOR_OPERATIONS_VECTOR_DOT_PRODUCT_HPP_

#include "../../meta/core_vector_meta.hpp"
#include "../concrete/core_vector_sharedmem_eigen.hpp"

namespace rompp{
namespace core{
namespace vec_ops{

template <typename vec_a_type,
	  typename vec_b_type,
	  typename res_t,
	  core::meta::enable_if_t<
	   details::traits<vec_a_type>::is_vector &&
	   details::traits<vec_a_type>::isEigen && 
	   details::traits<vec_b_type>::is_vector &&
	   details::traits<vec_b_type>::isEigen &&
	   std::is_same<
	     typename details::traits<vec_a_type>::scalar_t,
	     typename details::traits<vec_b_type>::scalar_t
	     >::value &&
	   std::is_same<
	     typename details::traits<vec_a_type>::scalar_t,
	     res_t>::value
	   > * = nullptr
	  >
void dot(const vec_a_type & vecA,
	 const vec_b_type & vecB,
	 res_t & result)
{
  assert(vecA.size() == vecB.size());
  result = vecA.data()->dot(*vecB.data());
}


} // end namespace vec_ops
} // end namespace core
}//end namespace rompp
#endif
