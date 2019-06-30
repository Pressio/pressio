
#ifndef ODE_IMPLICIT_POLICIES_FIND_IF_LEGITIMATE_IMPLICIT_RESIDUAL_POLICY_HPP_
#define ODE_IMPLICIT_POLICIES_FIND_IF_LEGITIMATE_IMPLICIT_RESIDUAL_POLICY_HPP_

#include "./ode_is_legitimate_implicit_residual_policy.hpp"

namespace rompp{ namespace ode{ namespace meta {

template<
  ImplicitEnum name,
  int nstates,
  typename state_t,
  typename residual_t,
  typename model_t,
  typename scalar_t,
  class ... Args2
  >
struct find_if_legitimate_implicit_residual_policy;

template<
  ImplicitEnum name,
  int nstates,
  typename state_t,
  typename residual_t,
  typename model_t,
  typename scalar_t
  >
struct find_if_legitimate_implicit_residual_policy<
  name, nstates, state_t, residual_t, model_t, scalar_t
  > : std::integral_constant<std::size_t, 0>{};


template<
  ImplicitEnum name,
  int nstates,
  typename state_t,
  typename residual_t,
  typename model_t,
  typename scalar_t,
  class Head, class ... Tail
  >
struct find_if_legitimate_implicit_residual_policy<
  name, nstates, state_t, residual_t, model_t, scalar_t,
  Head, Tail...
  >
  : std::conditional <
  is_legitimate_implicit_residual_policy<
    Head, name, nstates, state_t, residual_t,
    model_t, scalar_t
    >::type::value,
  std::integral_constant<std::size_t, 0>,
  std::integral_constant <
    std::size_t, 1 +
    find_if_legitimate_implicit_residual_policy
    <name, nstates, state_t, residual_t, model_t, scalar_t,
    Tail...>::type::value
    >
  >::type
{};


template <
  ImplicitEnum name,
  int nstates,
  typename state_t,
  typename residual_t,
  typename model_t,
  typename scalar_t,
  class... Args
  >
using find_if_legitimate_implicit_residual_policy_t =
  typename find_if_legitimate_implicit_residual_policy
  <name, nstates, state_t, residual_t, model_t, scalar_t,
   Args...>::type;


}}} // namespace rompp::ode::meta
#endif
