
#ifndef CORE_VECTOR_VECTOR_SUBTRACT_OPERATOR_SHAREDMEM_HPP_
#define CORE_VECTOR_VECTOR_SUBTRACT_OPERATOR_SHAREDMEM_HPP_

#include "core_vector_traits.hpp"
#include "core_vector_meta.hpp"
#include "core_vector_distributed_binary_expression_templates.hpp"
#include "core_vector_sharedmem_binary_expression_templates.hpp"


namespace rompp{ namespace core{
  
  
// T1: expre, T2: vector:
// example: a*3 - b
template <typename T1,
	  typename T2,
	  core::meta::enable_if_t<
  exprtemplates::is_admissible_vec_for_sharedmem_expression<T2>::value
	    > * = nullptr>
auto operator-(const T1 & u, const T2 & v) {
  using sc_t = typename core::details::traits<T2>::scalar_t;
  using ord_t = typename core::details::traits<T2>::ordinal_t;

  return core::exprtemplates::SharedMemVectorBinaryExp<
    core::exprtemplates::subtract_,T1, T2, sc_t, ord_t>(u, v);
}

//-----------------------------------------------------  
// T1: expre, T2: expr:
// example: a*3 - b*21
template <typename T1,
	  typename T2,
	  core::meta::enable_if_t<
  exprtemplates::is_sharedmem_vector_expression<T1>::value &&
  exprtemplates::is_sharedmem_vector_expression<T2>::value
	    > * = nullptr>
auto operator-(const T1 & u, const T2 & v) {
  using sc_t = typename T2::sc_type;
  using ord_t = typename T2::ord_type;

  return core::exprtemplates::SharedMemVectorBinaryExp<
    core::exprtemplates::subtract_,T1, T2, sc_t, ord_t>(u, v);
}
  
//-----------------------------------------------------
// T1: vector, T2: expr:
// example: a - b*21
template <typename T1,
	  typename T2,
	  core::meta::enable_if_t<
exprtemplates::is_admissible_vec_for_sharedmem_expression<T1>::value &&
exprtemplates::is_sharedmem_vector_expression<T2>::value
	    > * = nullptr>
auto operator-(const T1 & u, const T2 & v) {
  using sc_t = typename T2::sc_type;
  using ord_t = typename T2::ord_type;

  return core::exprtemplates::SharedMemVectorBinaryExp<
    core::exprtemplates::subtract_,T1, T2, sc_t, ord_t>(u, v);
}

  
}}//end namespace rompp::core
#endif
