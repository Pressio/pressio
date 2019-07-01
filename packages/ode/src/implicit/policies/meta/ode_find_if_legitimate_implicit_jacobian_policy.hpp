
#ifndef ODE_POLICIES_META_FIND_IF_LEGITIMATE_IMPLICIT_JACOBIAN_POLICY_HPP_
#define ODE_POLICIES_META_FIND_IF_LEGITIMATE_IMPLICIT_JACOBIAN_POLICY_HPP_

#include "./ode_is_legitimate_implicit_jacobian_policy.hpp"

namespace rompp{ namespace ode{ namespace meta {

template<
  ImplicitEnum name,
  typename state_t,
  typename jacobian_t,
  typename system_t,
  typename scalar_t,
  class ... Args2
  >
struct find_if_legitimate_implicit_jacobian_policy;

template<
  ImplicitEnum name,
  typename state_t,
  typename jacobian_t,
  typename system_t,
  typename scalar_t
  >
struct find_if_legitimate_implicit_jacobian_policy<
  name, state_t, jacobian_t, system_t, scalar_t
  > : std::integral_constant<std::size_t, 0>{};


template<
  ImplicitEnum name,
  typename state_t,
  typename jacobian_t,
  typename system_t,
  typename scalar_t,
  class Head, class ... Tail
  >
struct find_if_legitimate_implicit_jacobian_policy<
  name, state_t, jacobian_t, system_t, scalar_t,
  Head, Tail...
  >
  : std::conditional <
  is_legitimate_implicit_jacobian_policy
	  <Head, name,
	   state_t, jacobian_t,
	   system_t, scalar_t
           >::type::value,
  std::integral_constant<std::size_t, 0>,
  std::integral_constant <
    std::size_t, 1 +
    find_if_legitimate_implicit_jacobian_policy
    <name, state_t, jacobian_t, system_t, scalar_t, Tail...>::type::value
    >
  >::type
{};

template <ImplicitEnum name,
	  typename state_t,
	  typename jacobian_t,
	  typename system_t,
	  typename scalar_t,
	  class... Args>
using find_if_legitimate_implicit_jacobian_policy_t =
  typename find_if_legitimate_implicit_jacobian_policy<
  name, state_t, jacobian_t, system_t, scalar_t,
  Args...>::type;


}}} // namespace rompp::ode::meta
#endif
