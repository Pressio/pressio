
#ifndef ODE_POLICIES_META_IS_LEGITIMATE_IMPLICIT_JACOBIAN_POLICIES_HPP_
#define ODE_POLICIES_META_IS_LEGITIMATE_IMPLICIT_JACOBIAN_POLICIES_HPP_

#include "../base/ode_jacobian_policy_base.hpp"
#include "ode_has_jacobian_method_callable_wtih_four_args.hpp"
#include "ode_has_jacobian_method_callable_wtih_five_args.hpp"

namespace rompp{ namespace ode{ namespace meta {

template<
  typename T,
  ImplicitEnum name,
  typename state_t, typename jacobian_t,
  typename model_t, typename scalar_t,
  typename enable = void
  >
struct is_legitimate_implicit_jacobian_policy : std::false_type
{};


template<
  typename T,
  ImplicitEnum name,
  typename state_t, typename jacobian_t,
  typename model_t, typename scalar_t
  >
struct is_legitimate_implicit_jacobian_policy
<T, name, state_t, jacobian_t, model_t, scalar_t,
 ::rompp::mpl::enable_if_t<
   has_jacobian_callable_with_four_args<
     T, name, state_t, model_t, scalar_t>::value
   and
   has_jacobian_callable_with_five_args<
     T, name, state_t, jacobian_t,
     model_t, scalar_t>::value
   >
 > : std::true_type{};
//------------------------------------------------------------------

template<typename T, typename ... args>
using is_legitimate_implicit_euler_jacobian_policy =
  is_legitimate_implicit_jacobian_policy<
  T, ::rompp::ode::ImplicitEnum::Euler, args...>;

template<typename T, typename ... args>
using is_legitimate_implicit_bdf2_jacobian_policy =
  is_legitimate_implicit_jacobian_policy<
  T, ::rompp::ode::ImplicitEnum::BDF2, args...>;
//------------------------------------------------------------------


template<
  ImplicitEnum name,
  typename state_t, typename jacobian_t,
  typename model_t, typename scalar_t,
  class ... Args2
  >
struct find_legitimate_implicit_jacobian_policy;

template<
  ImplicitEnum name,
  typename state_t, typename jacobian_t,
  typename model_t, typename scalar_t
  >
struct find_legitimate_implicit_jacobian_policy<
  name, state_t, jacobian_t, model_t, scalar_t
  > : std::integral_constant<std::size_t, 0>{};


template<
  ImplicitEnum name,
  typename state_t, typename jacobian_t,
  typename model_t, typename scalar_t,
  class Head, class ... Tail
  >
struct find_legitimate_implicit_jacobian_policy<
  name, state_t, jacobian_t, model_t, scalar_t,
  Head, Tail...
  >
  : std::conditional <
  is_legitimate_implicit_jacobian_policy<Head, name,
					 state_t, jacobian_t,
					 model_t, scalar_t>::type::value,
  std::integral_constant<std::size_t, 0>,
  std::integral_constant <
    std::size_t, 1 +
    find_legitimate_implicit_jacobian_policy
    <name, state_t, jacobian_t, model_t, scalar_t, Tail...>::type::value
    >
  >::type
{};

template <ImplicitEnum name,
	  typename state_t, typename jacobian_t,
	  typename model_t, typename scalar_t,
	  class... Args>
using find_legitimate_implicit_jacobian_policy_t =
  typename find_legitimate_implicit_jacobian_policy<
  name, state_t, jacobian_t, model_t, scalar_t, Args...
  >::type;


}}} // namespace rompp::ode::meta
#endif
