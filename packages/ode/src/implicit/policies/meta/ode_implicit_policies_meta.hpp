
#ifndef ODE_POLICIES_META_IMPLICIT_POLICIES_META_HPP_
#define ODE_POLICIES_META_IMPLICIT_POLICIES_META_HPP_

#include "../base/ode_implicit_residual_policy_base.hpp"
#include "../base/ode_jacobian_policy_base.hpp"

namespace rompp{
namespace ode{
namespace meta {

// //-----------------------------------------------------------------
// // METAF FOR ADMISSIBLE JACOBIAN POLICY
// //-----------------------------------------------------------------
// template<typename policy_t, typename enable = void>
// struct is_legitimate_implicit_jacobian_policy : std::false_type{};
  
// template <typename policy_t>
// struct is_legitimate_implicit_jacobian_policy<
//   policy_t,
//   typename std::enable_if<
//     core::meta::publicly_inherits_from<
//       policy_t,
//       ode::policy::JacobianPolicyBase<policy_t>
//       >::value
//     >::type
//   > : std::true_type{};
// //-----------------------------------------------------------------


//-----------------------------------------------------------------
// ADMISSIBLE IMPLICIT JACOBIAN
//-----------------------------------------------------------------
template<::rompp::ode::ImplicitSteppersEnum whichone,
	  typename policy_t, typename enable = void>
struct is_legitimate_implicit_jacobian_policy : std::false_type{};
  
template <typename policy_t>
struct is_legitimate_implicit_jacobian_policy<
  ::rompp::ode::ImplicitSteppersEnum::Euler,
  policy_t,
  typename std::enable_if<
    core::meta::publicly_inherits_from<
      policy_t,
      ode::policy::JacobianPolicyBase<policy_t>
      >::value 
    >::type > : std::true_type{};

template <typename policy_t>
struct is_legitimate_implicit_jacobian_policy<
  ::rompp::ode::ImplicitSteppersEnum::BDF2,
  policy_t,
  typename std::enable_if<
    core::meta::publicly_inherits_from<
      policy_t,
      ode::policy::JacobianPolicyBase<policy_t>
      >::value 
    >::type > : std::true_type{};

  
//-----------------------------------------------------------------
// ADMISSIBLE IMPLICIT RESIDUAL
//-----------------------------------------------------------------
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
      ode::policy::ImplicitResidualPolicyBase<policy_t, 1, 0>
      >::value 
    >::type
  > : std::true_type{};
  
  
//-----------------------------------------------------------------
// METAF TO CHECK RESIDUAL POLICY IS STANDARD
//-----------------------------------------------------------------  
template<::rompp::ode::ImplicitSteppersEnum whichone,
	  typename policy_t, typename enable = void>
struct is_implicit_residual_standard_policy : std::false_type{};
  

  template <template <typename...> class policy_t, typename... Args>
struct is_implicit_residual_standard_policy<
  ::rompp::ode::ImplicitSteppersEnum::Euler,
  policy_t<Args...>,
  typename std::enable_if<
    std::is_same<policy_t<Args...>,
		 ode::policy::ImplicitEulerResidualStandardPolicy<
		 Args...>
		 >::value
    >::type > : std::true_type{};

  
template <template <typename...> class policy_t, typename... Args>
struct is_implicit_residual_standard_policy<
  ::rompp::ode::ImplicitSteppersEnum::BDF2,
  policy_t<Args...>,
  typename std::enable_if<
    std::is_same<policy_t<Args...>,
     ode::policy::ImplicitBDF2ResidualStandardPolicy<Args...>
		 >::value
    >::type > : std::true_type{};


//-----------------------------------------------------------------
// METAF TO CHECK JACOBIAN POLICY IS STANDARD
//-----------------------------------------------------------------  
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
  
  
  
} // namespace meta
} // namespace core
}//end namespace rompp
#endif
