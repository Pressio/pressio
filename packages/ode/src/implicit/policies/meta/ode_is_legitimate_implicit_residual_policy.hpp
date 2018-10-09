
#ifndef ODE_POLICIES_META_IS_LEGITIMATE_IMPLICIT_RESIDUAL_POLICIES_HPP_
#define ODE_POLICIES_META_IS_LEGITIMATE_IMPLICIT_RESIDUAL_POLICIES_HPP_

#include "../base/ode_implicit_residual_policy_base.hpp"

namespace rompp{ namespace ode{ namespace meta {

    
template<::rompp::ode::ImplicitSteppersEnum whichone,
	  typename policy_t, typename enable = void>
struct is_legitimate_implicit_residual_policy : std::false_type{};
  
template <typename policy_t>
struct is_legitimate_implicit_residual_policy<
  ::rompp::ode::ImplicitSteppersEnum::Euler,
  policy_t,
  typename std::enable_if<
    core::meta::publicly_inherits_from<
      policy_t,
      ode::policy::ImplicitResidualPolicyBase<policy_t, 1, 0>
      >::value 
    >::type
  > : std::true_type{};

template <typename policy_t>
struct is_legitimate_implicit_residual_policy<
  ::rompp::ode::ImplicitSteppersEnum::BDF2,
  policy_t,
  typename std::enable_if<
    core::meta::publicly_inherits_from<
      policy_t,
      ode::policy::ImplicitResidualPolicyBase<policy_t, 2, 0>
      >::value 
    >::type
  > : std::true_type{};


template<typename T>
using is_legitimate_implicit_euler_residual_policy =
  typename is_legitimate_implicit_residual_policy<
  ::rompp::ode::ImplicitSteppersEnum::Euler,T>::type;

template<typename T>
using is_legitimate_implicit_bdf2_residual_policy =
  typename is_legitimate_implicit_residual_policy<
  ::rompp::ode::ImplicitSteppersEnum::BDF2,T>::type;
  
  
  
}}} // namespace rompp::ode::meta
#endif
