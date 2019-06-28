
#ifndef ALGEBRA_VECTOR_VECTOR_SUBTRACT_OPERATOR_SHAREDMEM_HPP_
#define ALGEBRA_VECTOR_VECTOR_SUBTRACT_OPERATOR_SHAREDMEM_HPP_

#include "../algebra_vector_traits.hpp"
#include "../algebra_vector_meta.hpp"
#include "algebra_vector_sharedmem_binary_expression_templates.hpp"
#include "../../algebra_expression_templates_operators.hpp"

namespace rompp{ namespace algebra{

// T1: expre, T2: vector:
// example: a*3 - b
template <typename T1,
	  typename T2,
	  ::rompp::mpl::enable_if_t<
  meta::is_admissible_vec_for_sharedmem_expression<T2>::value
	    > * = nullptr>
auto operator-(const T1 & u, const T2 & v)
  -> decltype(
      algebra::exprtemplates::SharedMemVectorBinaryExp<
      algebra::exprtemplates::subtract_,
      T1, T2,
      typename algebra::details::traits<T2>::scalar_t,
      typename algebra::details::traits<T2>::ordinal_t>(u, v)
    )
{
  using sc_t = typename algebra::details::traits<T2>::scalar_t;
  using ord_t = typename algebra::details::traits<T2>::ordinal_t;

  return algebra::exprtemplates::SharedMemVectorBinaryExp<
    algebra::exprtemplates::subtract_, T1, T2, sc_t, ord_t>(u, v);
}

//-----------------------------------------------------
// T1: expre, T2: expr:
// example: a*3 - b*21
template <typename T1,
	  typename T2,
	  ::rompp::mpl::enable_if_t<
  exprtemplates::is_sharedmem_vector_expression<T1>::value &&
  exprtemplates::is_sharedmem_vector_expression<T2>::value
	    > * = nullptr>
auto operator-(const T1 & u, const T2 & v)
  -> decltype(
      algebra::exprtemplates::SharedMemVectorBinaryExp<
      algebra::exprtemplates::subtract_,
      T1, T2,
      typename T2::sc_type,
      typename T2::ord_type>(u, v)
    )
{
  using sc_t = typename T2::sc_type;
  using ord_t = typename T2::ord_type;

  return algebra::exprtemplates::SharedMemVectorBinaryExp<
    algebra::exprtemplates::subtract_, T1, T2, sc_t, ord_t>(u, v);
}

//-----------------------------------------------------
// T1: vector, T2: expr:
// example: a - b*21
template <typename T1,
	  typename T2,
	  ::rompp::mpl::enable_if_t<
meta::is_admissible_vec_for_sharedmem_expression<T1>::value &&
exprtemplates::is_sharedmem_vector_expression<T2>::value
	    > * = nullptr>
auto operator-(const T1 & u, const T2 & v)
  -> decltype(
      algebra::exprtemplates::SharedMemVectorBinaryExp<
      algebra::exprtemplates::subtract_,
      T1, T2,
      typename T2::sc_type,
      typename T2::ord_type>(u, v)
    )
{
  using sc_t = typename T2::sc_type;
  using ord_t = typename T2::ord_type;

  return algebra::exprtemplates::SharedMemVectorBinaryExp<
    algebra::exprtemplates::subtract_, T1, T2, sc_t, ord_t>(u, v);
}


}}//end namespace rompp::algebra
#endif
