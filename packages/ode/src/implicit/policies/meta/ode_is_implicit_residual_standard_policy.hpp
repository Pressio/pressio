
#ifndef ODE_IMPLICIT_POLICIES_IS_IMPLICIT_RESIDUAL_STD_POLICY_HPP_
#define ODE_IMPLICIT_POLICIES_IS_IMPLICIT_RESIDUAL_STD_POLICY_HPP_

#include "../base/ode_implicit_residual_policy_base.hpp"
#include "../standard/ode_implicit_residual_standard_policy.hpp"

namespace pressio{ namespace ode{ namespace meta {

template<
  ::pressio::ode::ImplicitEnum whichone,
  typename policy_t,
  typename enable = void
  >
struct is_implicit_residual_standard_policy : std::false_type{};

template <template <typename...> class policy_t, typename... Args>
struct is_implicit_residual_standard_policy<
  ::pressio::ode::ImplicitEnum::Euler,
  policy_t<Args...>,
  typename std::enable_if<
    std::is_same<
      policy_t<Args...>,
      ode::policy::ImplicitResidualStandardPolicy<Args...>
      >::value
    >::type > : std::true_type{};


template <template <typename...> class policy_t, typename... Args>
struct is_implicit_residual_standard_policy<
  ::pressio::ode::ImplicitEnum::BDF2,
  policy_t<Args...>,
  typename std::enable_if<
    std::is_same<
      policy_t<Args...>,
      ode::policy::ImplicitResidualStandardPolicy<Args...>
      >::value
    >::type > : std::true_type{};


template<template <typename...> class policy_t, typename... Args>
using is_implicit_euler_residual_standard_policy =
  typename is_implicit_residual_standard_policy<
  ::pressio::ode::ImplicitEnum::Euler, policy_t<Args...>>::type;

template<template <typename...> class policy_t, typename... Args>
using is_implicit_bdf2_residual_standard_policy =
  typename is_implicit_residual_standard_policy<
  ::pressio::ode::ImplicitEnum::BDF2, policy_t<Args...>>::type;



}}} // namespace pressio::ode::meta
#endif
