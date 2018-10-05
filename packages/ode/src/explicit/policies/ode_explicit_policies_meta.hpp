
#ifndef ODE_POLICIES_META_EXPLICIT_POLICIES_META_HPP_
#define ODE_POLICIES_META_EXPLICIT_POLICIES_META_HPP_

#include "ode_explicit_residual_standard_policy.hpp"

namespace rompp{
namespace ode{
namespace meta {

//-----------------------------------------------------------------
// METAF FOR ADMISSIBLE EXPLICIT RESIDUAL
//-----------------------------------------------------------------
template<typename policy_t, typename enable = void>
struct is_legitimate_explicit_residual_policy
  : std::false_type{};
  
template <typename policy_t>
struct is_legitimate_explicit_residual_policy<
  policy_t,
  typename std::enable_if<
    core::meta::publicly_inherits_from<
      policy_t,
      ode::policy::ExplicitResidualPolicyBase<policy_t>
      >::value
    >::type
  > : std::true_type{};

  
//-----------------------------------------------------------------
// CHECK IF RESIDUAL POLICY IS EULER STANDARD
//-----------------------------------------------------------------
template<typename policy_t, typename enable = void>
struct is_explicit_euler_residual_standard_policy
  : std::false_type{};

  
template <template <typename...> class policy_t,
	  typename ... Args>
struct is_explicit_euler_residual_standard_policy<
  policy_t<Args...>,
  typename std::enable_if<
    std::is_same<policy_t<Args...>,
		 ode::policy::ExplicitResidualStandardPolicy<
		   Args...>
		 >::value
    >::type
  > : std::true_type{};


//-----------------------------------------------------------------
// CHECK IF RESIDUAL POLICY IS RK4 STANDARD
//-----------------------------------------------------------------
template<typename policy_t, typename enable = void>
struct is_explicit_runge_kutta4_residual_standard_policy
  : std::false_type{};

template <template <typename...> class policy_t,
	  typename ... Args>
struct is_explicit_runge_kutta4_residual_standard_policy<
  policy_t<Args...>,
  typename std::enable_if<
    std::is_same<policy_t<Args...>,
		 ode::policy::ExplicitResidualStandardPolicy<
		   Args...>
		 >::value
    >::type
  > : std::true_type{};
//----------------------------------------------------------------

  
} // namespace meta
} // namespace core
}//end namespace rompp
#endif
