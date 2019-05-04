
#ifndef ODE_POLICIES_META_HAS_RESIDUAL_METHOD_CALLABLE_WITH_FIVE_ARGS_HPP_
#define ODE_POLICIES_META_HAS_RESIDUAL_METHOD_CALLABLE_WITH_FIVE_ARGS_HPP_

#include "../base/ode_implicit_residual_policy_base.hpp"

namespace rompp{ namespace ode{ namespace meta {

template <
  typename T,
  ::rompp::ode::ImplicitEnum odeMethod,
  int numAuxStates,
  typename state_t,
  typename model_t,
  typename scalar_t,
  typename enable = void
  >
struct has_residual_callable_with_five_args
  : std::false_type{};


template <
  typename T,
  ::rompp::ode::ImplicitEnum odeMethod,
  int numAuxStates,
  typename state_t,
  typename model_t,
  typename scalar_t
  >
struct has_residual_callable_with_five_args<
  T, odeMethod, numAuxStates, state_t, model_t, scalar_t,
  ::rompp::mpl::enable_if_t<
    !std::is_void<
      decltype(std::declval<T>().template operator()<
	   odeMethod,
	   numAuxStates>( std::declval<const state_t &>(),
			  std::declval<const std::array<state_t, numAuxStates> &>(),
			  std::declval<const model_t&>(),
			  std::declval<scalar_t>(),
			  std::declval<scalar_t>()
			  )
	       )
      >::value
    >
  > : std::true_type
{};


}}} // namespace rompp::ode::meta
#endif
