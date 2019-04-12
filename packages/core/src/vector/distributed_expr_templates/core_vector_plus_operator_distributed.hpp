
#ifndef CORE_VECTOR_VECTOR_PLUS_OPERATOR_DISTRIBUTED_HPP_
#define CORE_VECTOR_VECTOR_PLUS_OPERATOR_DISTRIBUTED_HPP_

#include "../core_vector_traits.hpp"
#include "../core_vector_meta.hpp"
#include "core_vector_distributed_binary_expression_templates.hpp"
#include "../../core_expression_templates_operators.hpp"

namespace rompp{ namespace core{

// T1: expre, T2: vector:
// example: a*3 + b
template <typename T1,
	  typename T2,
	  ::rompp::mpl::enable_if_t<
  meta::is_admissible_vec_for_dist_expression<T2>::value
	    > * = nullptr>
auto operator+(const T1 & u, const T2 & v)
  -> decltype(
      core::exprtemplates::DistributedVectorBinaryExp<
      core::exprtemplates::plus_,
      T1, T2,
      typename core::details::traits<T2>::scalar_t,
      typename core::details::traits<T2>::local_ordinal_t>(u, v)
    )
{
  using sc_t = typename core::details::traits<T2>::scalar_t;
  using LO_t = typename core::details::traits<T2>::local_ordinal_t;

  return core::exprtemplates::DistributedVectorBinaryExp<
    core::exprtemplates::plus_,T1, T2, sc_t, LO_t>(u, v);
}

//-----------------------------------------------------
// T1: expre, T2: expr:
// example: a*3 + b*21
template <typename T1,
	  typename T2,
	  ::rompp::mpl::enable_if_t<
  exprtemplates::is_distributed_vector_expression<T1>::value &&
  exprtemplates::is_distributed_vector_expression<T2>::value
	    > * = nullptr>
auto operator+(const T1 & u, const T2 & v)
  -> decltype(
      core::exprtemplates::DistributedVectorBinaryExp<
      core::exprtemplates::plus_,
      T1, T2,
      typename T2::sc_type,
      typename T2::LO_type>(u, v)
    )
{
  using sc_t = typename T2::sc_type;
  using LO_t = typename T2::LO_type;

  return core::exprtemplates::DistributedVectorBinaryExp<
    core::exprtemplates::plus_, T1, T2, sc_t, LO_t>(u, v);
}

//-----------------------------------------------------
// T1: vector, T2: expr:
// example: a + b*21
template <typename T1,
	  typename T2,
	  ::rompp::mpl::enable_if_t<
  meta::is_admissible_vec_for_dist_expression<T1>::value &&
  exprtemplates::is_distributed_vector_expression<T2>::value
	    > * = nullptr>
auto operator+(const T1 & u, const T2 & v)
  -> decltype(
      core::exprtemplates::DistributedVectorBinaryExp<
      core::exprtemplates::plus_,
      T1, T2,
      typename T2::sc_type,
      typename T2::LO_type>(u, v)
    )
{
  using sc_t = typename T2::sc_type;
  using LO_t = typename T2::LO_type;

  return core::exprtemplates::DistributedVectorBinaryExp<
    core::exprtemplates::plus_, T1, T2, sc_t, LO_t>(u, v);
}

}}//end namespace rompp::core
#endif
