
#ifndef ALGEBRA_VECTOR_VECTOR_TIMES_OPERATOR_SHAREDMEM_HPP_
#define ALGEBRA_VECTOR_VECTOR_TIMES_OPERATOR_SHAREDMEM_HPP_

#include "../algebra_vector_traits.hpp"
#include "../algebra_vector_meta.hpp"
#include "algebra_vector_sharedmem_binary_expression_templates.hpp"
#include "../../algebra_expression_templates_operators.hpp"

namespace rompp{ namespace algebra{

// T1: scalar, T2: vector:
// example: 3.*a
template <typename T1, typename T2,
	    ::rompp::mpl::enable_if_t<
  std::is_scalar<T1>::value &&
  meta::is_admissible_vec_for_sharedmem_expression<T2>::value
	      > * = nullptr>
auto operator*(T1 u, const T2 & v)
  -> decltype(
      algebra::exprtemplates::SharedMemVectorBinaryExp<
      algebra::exprtemplates::times_,
      T2, T1, T1, typename algebra::details::traits<T2>::ordinal_t>(v,u)
    )
{
  using sc_t = T1;
  using vec_sc_t = typename algebra::details::traits<T2>::scalar_t;
  static_assert(std::is_same<sc_t, vec_sc_t>::value, "");
  using ord_t = typename algebra::details::traits<T2>::ordinal_t;

  return algebra::exprtemplates::SharedMemVectorBinaryExp<
    algebra::exprtemplates::times_, T2, sc_t, sc_t, ord_t>(v,u);
}
//-----------------------------------------------------

// T1: vector, T2: scalar:
// example: a*3
template <typename T1, typename T2,
	    ::rompp::mpl::enable_if_t<
  std::is_scalar<T2>::value &&
  meta::is_admissible_vec_for_sharedmem_expression<T1>::value
	      > * = nullptr>
auto operator*(const T1 & u, T2 v)
  -> decltype(
      algebra::exprtemplates::SharedMemVectorBinaryExp<
      algebra::exprtemplates::times_,
      T1, T2, T2, typename algebra::details::traits<T1>::ordinal_t>(u,v)
    )
{
  using sc_t = T2;
  using vec_sc_t = typename algebra::details::traits<T1>::scalar_t;
  static_assert(std::is_same<sc_t, vec_sc_t>::value, "");
  using ord_t = typename algebra::details::traits<T1>::ordinal_t;

  return algebra::exprtemplates::SharedMemVectorBinaryExp<
    algebra::exprtemplates::times_, T1, sc_t, sc_t, ord_t>(u,v);
}
//-----------------------------------------------------

// T1: expre, T2: scalar:
// example: (a + b)*2
template <typename T1,
	  typename T2,
	  ::rompp::mpl::enable_if_t<
  exprtemplates::is_sharedmem_vector_expression<T1>::value &&
  std::is_scalar<T2>::value
	    > * = nullptr>
auto operator*(const T1 & u, T2 v)
  -> decltype(
      algebra::exprtemplates::SharedMemVectorBinaryExp<
      algebra::exprtemplates::times_,
      T1, typename T1::sc_type, typename T1::sc_type,
      typename T1::ord_type>(u,v)
    )
{
  using sc_t = typename T1::sc_type;
  using ord_t = typename T1::ord_type;

  return algebra::exprtemplates::SharedMemVectorBinaryExp<
    algebra::exprtemplates::times_, T1, sc_t, sc_t, ord_t>(u, v);
}
//-----------------------------------------------------

// T1: scalar, T2: expr:
// example: 2*(a + b)
template <typename T1,
	  typename T2,
	  ::rompp::mpl::enable_if_t<
  std::is_scalar<T1>::value &&
  exprtemplates::is_sharedmem_vector_expression<T2>::value
	    > * = nullptr>
auto operator*(T1 u, const T2 & v)
  -> decltype(
      algebra::exprtemplates::SharedMemVectorBinaryExp<
      algebra::exprtemplates::times_,
      T2, typename T2::sc_type, typename T2::sc_type,
      typename T2::ord_type>(v,u)
    )
{
  using sc_t = typename T2::sc_type;
  using ord_t = typename T2::ord_type;

  return algebra::exprtemplates::SharedMemVectorBinaryExp<
    algebra::exprtemplates::times_, T2, sc_t, sc_t, ord_t>(v,u);
}
//-----------------------------------------------------

// T1: vector, T2: vector
// example: u*v
template <typename T1, typename T2,
	    ::rompp::mpl::enable_if_t<
  meta::is_admissible_vec_for_sharedmem_expression<T1>::value and
  meta::is_admissible_vec_for_sharedmem_expression<T2>::value
	      > * = nullptr>
auto operator*(const T1 & u, const T2 & v)
  -> decltype(
      algebra::exprtemplates::SharedMemVectorBinaryExp<
      algebra::exprtemplates::times_,
      T1, T2,
      typename algebra::details::traits<T1>::scalar_t,
      typename algebra::details::traits<T1>::ordinal_t>(u,v)
    )
{
  using vec_u_sc_t = typename algebra::details::traits<T1>::scalar_t;
  using vec_v_sc_t = typename algebra::details::traits<T2>::scalar_t;
  static_assert(std::is_same<vec_u_sc_t, vec_v_sc_t>::value, "");
  using ord_t = typename algebra::details::traits<T1>::ordinal_t;

  return algebra::exprtemplates::SharedMemVectorBinaryExp<
    algebra::exprtemplates::times_, T1, T2, vec_u_sc_t, ord_t>(u,v);
}
//-----------------------------------------------------

// T1: expre, T2: vector:
// example: (a*b)*c
template <typename T1,
	  typename T2,
	  ::rompp::mpl::enable_if_t<
  exprtemplates::is_sharedmem_vector_expression<T1>::value &&
  meta::is_admissible_vec_for_sharedmem_expression<T2>::value
	    > * = nullptr>
auto operator*(const T1 & u, const T2 v)
  -> decltype(
      algebra::exprtemplates::SharedMemVectorBinaryExp<
      algebra::exprtemplates::times_,
      T1, T2,
      typename T1::sc_type, typename T1::ord_type>(u,v)
    )
{
  using sc_t = typename T1::sc_type;
  using ord_t = typename T1::ord_type;
  using vec_v_sc_t = typename algebra::details::traits<T2>::scalar_t;
  static_assert(std::is_same<sc_t, vec_v_sc_t>::value, "");

  return algebra::exprtemplates::SharedMemVectorBinaryExp<
    algebra::exprtemplates::times_, T1, T2, sc_t, ord_t>(u, v);
}
//-----------------------------------------------------


}}//end namespace rompp::algebra
#endif // ALGEBRA_VECTOR_VECTOR_TIMES_OPERATOR_SHAREDMEM_HPP_
