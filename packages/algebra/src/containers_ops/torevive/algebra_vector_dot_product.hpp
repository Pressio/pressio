
#ifndef ALGEBRA_VECTOR_OPERATIONS_VECTOR_DOT_PRODUCT_HPP_
#define ALGEBRA_VECTOR_OPERATIONS_VECTOR_DOT_PRODUCT_HPP_

#include "../../meta/algebra_vector_meta.hpp"
#include "../concrete/algebra_vector_sharedmem_eigen.hpp"

namespace rompp{
namespace algebra{
namespace vec_ops{

template <typename vec_a_type,
	  typename vec_b_type,
	  typename res_t,
	  ::rompp::mpl::enable_if_t<
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
} // end namespace algebra
}//end namespace rompp
#endif
