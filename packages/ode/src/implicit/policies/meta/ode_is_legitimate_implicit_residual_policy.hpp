
#ifndef ODE_IMPLICIT_POLICIES_IS_LEGITIMATE_IMPLICIT_RESIDUAL_POLICIES_HPP_
#define ODE_IMPLICIT_POLICIES_IS_LEGITIMATE_IMPLICIT_RESIDUAL_POLICIES_HPP_

#include "../base/ode_implicit_residual_policy_base.hpp"
#include "./ode_has_residual_method_callable_wtih_five_args.hpp"
#include "./ode_has_residual_method_callable_wtih_six_args.hpp"

namespace rompp{ namespace ode{ namespace meta {

template<
  typename T,
  ImplicitEnum name, int nstates,
  typename state_t, typename residual_t,
  typename model_t, typename scalar_t,
  typename enable = void
  >
struct is_legitimate_implicit_residual_policy : std::false_type{};

template<
  typename T,
  ImplicitEnum name, int nstates,
  typename state_t, typename residual_t,
  typename model_t, typename scalar_t
  >
struct is_legitimate_implicit_residual_policy<
  T, name, nstates, state_t, residual_t, model_t, scalar_t,
  ::rompp::mpl::enable_if_t<
    has_residual_callable_with_six_args<
      T, name, nstates, state_t, residual_t, model_t, scalar_t
      >::value and
    has_residual_callable_with_five_args<
      T, name, nstates, state_t, model_t, scalar_t
      >::value
    >
  > : std::true_type{};
//------------------------------------------------------------------


template<typename T, typename ... args>
using is_legitimate_implicit_euler_residual_policy =
  is_legitimate_implicit_residual_policy<
  T, ImplicitEnum::Euler, 1, args...>;

template<typename T, typename ... args>
using is_legitimate_implicit_bdf2_residual_policy =
  is_legitimate_implicit_residual_policy<
  T, ImplicitEnum::BDF2, 2, args...>;
//------------------------------------------------------------------


template<
  ImplicitEnum name, int nstates,
  typename state_t, typename residual_t,
  typename model_t, typename scalar_t,
  class ... Args2
  >
struct find_legitimate_implicit_residual_policy;

template<
  ImplicitEnum name, int nstates,
  typename state_t, typename residual_t,
  typename model_t, typename scalar_t
  >
struct find_legitimate_implicit_residual_policy<
  name, nstates, state_t, residual_t, model_t, scalar_t
  > : std::integral_constant<std::size_t, 0>{};


template<
  ImplicitEnum name, int nstates,
  typename state_t, typename residual_t,
  typename model_t, typename scalar_t,
  class Head, class ... Tail
  >
struct find_legitimate_implicit_residual_policy<
  name, nstates, state_t, residual_t, model_t, scalar_t,
  Head, Tail...
  >
  : std::conditional <
  is_legitimate_implicit_residual_policy<Head, name, nstates,
					 state_t, residual_t,
					 model_t, scalar_t>::type::value,
  std::integral_constant<std::size_t, 0>,
  std::integral_constant <
    std::size_t, 1 +
    find_legitimate_implicit_residual_policy
    <name, nstates, state_t, residual_t, model_t, scalar_t, Tail...>::type::value
    >
  >::type
{};

template <ImplicitEnum name, int nstates,
	  typename state_t, typename residual_t,
	  typename model_t, typename scalar_t,
	  class... Args>
using find_legitimate_implicit_residual_policy_t =
  typename find_legitimate_implicit_residual_policy
  <name, nstates, state_t, residual_t, model_t, scalar_t, Args...>::type;


}}} // namespace rompp::ode::meta
#endif
