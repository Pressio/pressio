
#ifndef ODE_EXPLICIT_POLICIES_IS_EXPLICIT_RK4_VELOCITY_STANDARD_POLICY_HPP_
#define ODE_EXPLICIT_POLICIES_IS_EXPLICIT_RK4_VELOCITY_STANDARD_POLICY_HPP_

#include "ode_explicit_velocity_standard_policy.hpp"

namespace rompp{ namespace ode{ namespace meta {

template<typename policy_t, typename enable = void>
struct is_explicit_runge_kutta4_velocity_standard_policy
  : std::false_type{};

template <template <typename...> class policy_t,
	  typename ... Args>
struct is_explicit_runge_kutta4_velocity_standard_policy<
  policy_t<Args...>,
  typename std::enable_if<
    std::is_same<policy_t<Args...>,
		 ode::policy::ExplicitVelocityStandardPolicy<
		   Args...>
		 >::value
    >::type
  > : std::true_type{};

}}}//end namespace rompp::containers::meta
#endif
