
#ifndef ODE_POLICIES_META_IS_IMPLICIT_RESIDUAL_STD_POLICY_HPP_
#define ODE_POLICIES_META_IS_IMPLICIT_RESIDUAL_STD_POLICY_HPP_

#include "../base/ode_jacobian_policy_base.hpp"
#include "../standard/ode_implicit_euler_jacobian_standard_policy.hpp"
#include "../standard/ode_implicit_bdf2_jacobian_standard_policy.hpp"

namespace rompp{ namespace ode{ namespace meta {

      
template<::rompp::ode::ImplicitSteppersEnum whichone,
	  typename policy_t, typename enable = void>
struct is_implicit_jacobian_standard_policy : std::false_type{};
  
  template <template <typename...> class policy_t, typename... Args>
struct is_implicit_jacobian_standard_policy<
  ::rompp::ode::ImplicitSteppersEnum::Euler,
  policy_t<Args...>,
  typename std::enable_if<
    std::is_same<policy_t<Args...>,
		 ode::policy::ImplicitEulerJacobianStandardPolicy<
		 Args...>
		 >::value
    >::type > : std::true_type{};

  
template <template <typename...> class policy_t, typename... Args>
struct is_implicit_jacobian_standard_policy<
  ::rompp::ode::ImplicitSteppersEnum::BDF2,
  policy_t<Args...>,
  typename std::enable_if<
    std::is_same<policy_t<Args...>,
		 ode::policy::ImplicitBDF2JacobianStandardPolicy<
		 Args...>
		 >::value
    >::type > : std::true_type{};


template<template <typename...> class policy_t, typename... Args>
using is_implicit_euler_jacobian_standard_policy =
  typename is_implicit_jacobian_standard_policy<
  ::rompp::ode::ImplicitSteppersEnum::Euler,policy_t<Args...>>::type;

template<template <typename...> class policy_t, typename... Args>
using is_implicit_bdf2_jacobian_standard_policy =
  typename is_implicit_jacobian_standard_policy<
  ::rompp::ode::ImplicitSteppersEnum::BDF2,policy_t<Args...>>::type;
  
  
}}} // namespace rompp::ode::meta
#endif
