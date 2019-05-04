
#ifndef ODE_IMPLICIT_POLICIES_META_HAS_JACOBIAN_CALLABLE_WITH_FIVE_ARGS_HPP_
#define ODE_IMPLICIT_POLICIES_META_HAS_JACOBIAN_CALLABLE_WITH_FIVE_ARGS_HPP_

#include "../base/ode_jacobian_policy_base.hpp"

namespace rompp{ namespace ode{ namespace meta {

template <
  typename T,
  ::rompp::ode::ImplicitEnum odeMethod,
  typename state_t, typename jacobian_t, typename model_t, typename scalar_t,
  typename enable = void
  >
struct has_jacobian_callable_with_five_args : std::false_type
{};


template <
  typename T, ::rompp::ode::ImplicitEnum odeMethod,
  typename state_t, typename jacobian_t, typename model_t, typename scalar_t
  >
struct has_jacobian_callable_with_five_args<
  T, odeMethod, state_t, jacobian_t, model_t, scalar_t,
  ::rompp::mpl::enable_if_t<
    std::is_void<
      decltype(std::declval<T>().template operator()<
	       odeMethod>( std::declval<const state_t &>(),
			   std::declval<jacobian_t &>(),
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
