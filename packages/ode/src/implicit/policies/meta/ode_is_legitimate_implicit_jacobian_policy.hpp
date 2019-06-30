
#ifndef ODE_POLICIES_META_IS_LEGITIMATE_IMPLICIT_JACOBIAN_POLICIES_HPP_
#define ODE_POLICIES_META_IS_LEGITIMATE_IMPLICIT_JACOBIAN_POLICIES_HPP_

#include "../base/ode_jacobian_policy_base.hpp"

namespace rompp{ namespace ode{ namespace meta {

template<
  typename T,
  ImplicitEnum name,
  typename state_t,
  typename jacobian_t,
  typename model_t,
  typename scalar_t,
  typename enable = void
  >
struct is_legitimate_implicit_jacobian_policy : std::false_type
{};


template<
  typename T,
  ImplicitEnum name,
  typename state_t,
  typename jacobian_t,
  typename model_t,
  typename scalar_t
  >
struct is_legitimate_implicit_jacobian_policy
<T, name, state_t, jacobian_t, model_t, scalar_t,
 ::rompp::mpl::enable_if_t<
   std::is_same<
     jacobian_t,
     decltype
     (
      std::declval<T>().template operator()
      <
      name
      >( std::declval<const state_t &>(),
	 std::declval<const model_t&>(),
	 std::declval<scalar_t>(),
	 std::declval<scalar_t>()
	 )
      )
     >::value
   and
   std::is_void<
     decltype
     (
      std::declval<T>().template operator()
      <
      name
      >( std::declval<const state_t &>(),
	 std::declval<jacobian_t &>(),
	 std::declval<const model_t&>(),
	 std::declval<scalar_t>(),
	 std::declval<scalar_t>()
	 )
      )
   >::value
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

}}} // namespace rompp::ode::meta
#endif
