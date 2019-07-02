
#ifndef CONTAINERS_VECTOR_VECTOR_TIMES_OPERATOR_DISTRIBUTED_HPP_
#define CONTAINERS_VECTOR_VECTOR_TIMES_OPERATOR_DISTRIBUTED_HPP_

#include "../containers_vector_traits.hpp"
#include "../containers_vector_meta.hpp"
#include "containers_vector_distributed_binary_expression_templates.hpp"
#include "../../containers_expression_templates_operators.hpp"

namespace pressio{ namespace containers{

// T1: scalar, T2: vector:
// example: 3.*a
template <typename T1, typename T2,
	    ::pressio::mpl::enable_if_t<
  std::is_scalar<T1>::value &&
  meta::is_admissible_vec_for_dist_expression<T2>::value
	      > * = nullptr>
auto operator*(T1 u, const T2 & v)
  -> decltype(
      containers::exprtemplates::DistributedVectorBinaryExp<
      containers::exprtemplates::times_,
      T2, T1, T1, typename containers::details::traits<T2>::local_ordinal_t>(v, u)
    )
{
  using sc_t = T1;
  using vec_sc_t = typename containers::details::traits<T2>::scalar_t;
  static_assert(std::is_same<sc_t, vec_sc_t>::value, "");
  using LO_t = typename containers::details::traits<T2>::local_ordinal_t;

  return containers::exprtemplates::DistributedVectorBinaryExp<
    containers::exprtemplates::times_, T2, sc_t, sc_t, LO_t>(v,u);
}
//-----------------------------------------------------


// T1: vector, T2: scalar:
// example: a*3
template <typename T1, typename T2,
	    ::pressio::mpl::enable_if_t<
  std::is_scalar<T2>::value &&
  meta::is_admissible_vec_for_dist_expression<T1>::value
	      > * = nullptr>
auto operator*(const T1 & u, T2 v)
  -> decltype(
      containers::exprtemplates::DistributedVectorBinaryExp<
      containers::exprtemplates::times_,
      T1, T2, T2,
      typename containers::details::traits<T1>::local_ordinal_t>(u, v)
    )
{
  using sc_t = T2;
  using vec_sc_t = typename containers::details::traits<T1>::scalar_t;
  static_assert(std::is_same<sc_t, vec_sc_t>::value, "");
  using LO_t = typename containers::details::traits<T1>::local_ordinal_t;

  return containers::exprtemplates::DistributedVectorBinaryExp<
    containers::exprtemplates::times_, T1, sc_t, sc_t, LO_t>(u,v);
}
//-----------------------------------------------------


// T1: expre, T2: scalar:
// example: (a + b)*2
template <typename T1,
	  typename T2,
	  ::pressio::mpl::enable_if_t<
  exprtemplates::is_distributed_vector_expression<T1>::value &&
  std::is_scalar<T2>::value
	    > * = nullptr>
auto operator*(const T1 & u, T2 v)
  -> decltype(
      containers::exprtemplates::DistributedVectorBinaryExp<
      containers::exprtemplates::times_,
      T1, typename T1::sc_type, typename T1::sc_type,
      typename T1::LO_type>(u, v)
    )
{
  using sc_t = typename T1::sc_type;
  using LO_t = typename T1::LO_type;
  return containers::exprtemplates::DistributedVectorBinaryExp<
    containers::exprtemplates::times_, T1, sc_t, sc_t, LO_t>(u, v);
}
//-----------------------------------------------------


// T1: scalar, T2: expr:
// example: 2*(a + b)
template <typename T1,
	  typename T2,
	  ::pressio::mpl::enable_if_t<
  std::is_scalar<T1>::value &&
  exprtemplates::is_distributed_vector_expression<T2>::value
	    > * = nullptr>
auto operator*(T1 u, const T2 & v)
  -> decltype(
      containers::exprtemplates::DistributedVectorBinaryExp<
      containers::exprtemplates::times_,
      T2, typename T2::sc_type, typename T2::sc_type,
      typename T2::LO_type>(v,u)
    )
{
  using sc_t = typename T2::sc_type;
  using LO_t = typename T2::LO_type;
  return containers::exprtemplates::DistributedVectorBinaryExp<
    containers::exprtemplates::times_, T2, sc_t, sc_t, LO_t>(v,u);
}


}}//end namespace pressio::containers

#endif // CONTAINERS_VECTOR_VECTOR_TIMES_OPERATOR_DISTRIBUTED_HPP_
