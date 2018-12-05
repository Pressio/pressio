
#ifndef ODE_POLICIES_META_IS_LEGITIMATE_IMPLICIT_JACOBIAN_POLICIES_HPP_
#define ODE_POLICIES_META_IS_LEGITIMATE_IMPLICIT_JACOBIAN_POLICIES_HPP_

#include "../base/ode_jacobian_policy_base.hpp"

namespace rompp{ namespace ode{ namespace meta {

template<::rompp::ode::ImplicitEnum whichone,
	  typename policy_t, typename enable = void>
struct is_legitimate_implicit_jacobian_policy : std::false_type{};
  
template <typename policy_t>
struct is_legitimate_implicit_jacobian_policy<
  ::rompp::ode::ImplicitEnum::Euler,
  policy_t,
  typename std::enable_if<
    core::meta::publicly_inherits_from<
      policy_t,
      ode::policy::JacobianPolicyBase<policy_t>
      >::value 
    >::type > : std::true_type{};

template <typename policy_t>
struct is_legitimate_implicit_jacobian_policy<
  ::rompp::ode::ImplicitEnum::BDF2,
  policy_t,
  typename std::enable_if<
    core::meta::publicly_inherits_from<
      policy_t,
      ode::policy::JacobianPolicyBase<policy_t>
      >::value 
    >::type > : std::true_type{};


template<typename T>
using is_legitimate_implicit_euler_jacobian_policy =
  is_legitimate_implicit_jacobian_policy<
  ::rompp::ode::ImplicitEnum::Euler,T>;

template<typename T>
using is_legitimate_implicit_bdf2_jacobian_policy =
  is_legitimate_implicit_jacobian_policy<
  ::rompp::ode::ImplicitEnum::BDF2,T>;
  
   
}}} // namespace rompp::ode::meta
#endif
