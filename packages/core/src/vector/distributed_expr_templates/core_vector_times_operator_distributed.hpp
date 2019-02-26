
#ifndef CORE_VECTOR_VECTOR_TIMES_OPERATOR_DISTRIBUTED_HPP_
#define CORE_VECTOR_VECTOR_TIMES_OPERATOR_DISTRIBUTED_HPP_

#include "../core_vector_traits.hpp"
#include "../core_vector_meta.hpp"
#include "core_vector_distributed_binary_expression_templates.hpp"
#include "../../core_expression_templates_operators.hpp"

namespace rompp{ namespace core{

// T1: scalar, T2: vector:
// example: 3.*a
template <typename T1, typename T2,
	    core::meta::enable_if_t<
  std::is_scalar<T1>::value &&
  meta::is_admissible_vec_for_dist_expression<T2>::value
	      > * = nullptr>
auto operator*(T1 u, const T2 & v)
  -> decltype(
      core::exprtemplates::DistributedVectorBinaryExp<
      core::exprtemplates::times_,
      T2, T1, T1, typename core::details::traits<T2>::local_ordinal_t>(v, u)
    )
{
  using sc_t = T1;
  using vec_sc_t = typename core::details::traits<T2>::scalar_t;
  static_assert(std::is_same<sc_t, vec_sc_t>::value, "");
  using LO_t = typename core::details::traits<T2>::local_ordinal_t;

  return core::exprtemplates::DistributedVectorBinaryExp<
    core::exprtemplates::times_, T2, sc_t, sc_t, LO_t>(v,u);
}
//-----------------------------------------------------


// T1: vector, T2: scalar:
// example: a*3
template <typename T1, typename T2,
	    core::meta::enable_if_t<
  std::is_scalar<T2>::value &&
  meta::is_admissible_vec_for_dist_expression<T1>::value
	      > * = nullptr>
auto operator*(const T1 & u, T2 v)
  -> decltype(
      core::exprtemplates::DistributedVectorBinaryExp<
      core::exprtemplates::times_,
      T1, T2, T2,
      typename core::details::traits<T1>::local_ordinal_t>(u, v)
    )
{
  using sc_t = T2;
  using vec_sc_t = typename core::details::traits<T1>::scalar_t;
  static_assert(std::is_same<sc_t, vec_sc_t>::value, "");
  using LO_t = typename core::details::traits<T1>::local_ordinal_t;

  return core::exprtemplates::DistributedVectorBinaryExp<
    core::exprtemplates::times_, T1, sc_t, sc_t, LO_t>(u,v);
}
//-----------------------------------------------------


// T1: expre, T2: scalar:
// example: (a + b)*2
template <typename T1,
	  typename T2,
	  core::meta::enable_if_t<
  exprtemplates::is_distributed_vector_expression<T1>::value &&
  std::is_scalar<T2>::value
	    > * = nullptr>
auto operator*(const T1 & u, T2 v)
  -> decltype(
      core::exprtemplates::DistributedVectorBinaryExp<
      core::exprtemplates::times_,
      T1, typename T1::sc_type, typename T1::sc_type,
      typename T1::LO_type>(u, v)
    )
{
  using sc_t = typename T1::sc_type;
  using LO_t = typename T1::LO_type;
  return core::exprtemplates::DistributedVectorBinaryExp<
    core::exprtemplates::times_, T1, sc_t, sc_t, LO_t>(u, v);
}
//-----------------------------------------------------


// T1: scalar, T2: expr:
// example: 2*(a + b)
template <typename T1,
	  typename T2,
	  core::meta::enable_if_t<
  std::is_scalar<T1>::value &&
  exprtemplates::is_distributed_vector_expression<T2>::value
	    > * = nullptr>
auto operator*(T1 u, const T2 & v)
  -> decltype(
      core::exprtemplates::DistributedVectorBinaryExp<
      core::exprtemplates::times_,
      T2, typename T2::sc_type, typename T2::sc_type,
      typename T2::LO_type>(v,u)
    )
{
  using sc_t = typename T2::sc_type;
  using LO_t = typename T2::LO_type;
  return core::exprtemplates::DistributedVectorBinaryExp<
    core::exprtemplates::times_, T2, sc_t, sc_t, LO_t>(v,u);
}


}}//end namespace rompp::core

#endif // CORE_VECTOR_VECTOR_TIMES_OPERATOR_DISTRIBUTED_HPP_
