
#ifndef ODE_EXPLICIT_RUNGE_KUTTA4_POLICIES_META_HPP_
#define ODE_EXPLICIT_RUNGE_KUTTA4_POLICIES_META_HPP_

#include "ode_explicit_policies_meta.hpp"
#include "../standard/ode_explicit_residual_standard_policy.hpp"

namespace ode{
namespace meta {

//-----------------------------------------------------------------
// METAF TO CHECK RESIDUAL POLICY IS STANDARD
//-----------------------------------------------------------------
template<typename policy_t, typename enable = void>
struct isExplicitRungeKutta4ResidualStandardPolicy
  : std::false_type{};

template <template <typename...> class policy_t,
	  typename state_type,
	  typename residual_type,
	  typename model_type,
	  typename time_type,
	  typename sizer_type>
struct isExplicitRungeKutta4ResidualStandardPolicy<
  policy_t<state_type, residual_type,
	   model_type, time_type, sizer_type>,
  typename std::enable_if<
    std::is_same<policy_t<state_type, residual_type,
			  model_type, time_type, sizer_type>,
		 ode::policy::explicitResidualStandardPolicy<
		   state_type, residual_type,
		   model_type, time_type, sizer_type>
		 >::value
    >::type
  > : std::true_type{};
//----------------------------------------------------------------

} // namespace meta
} // namespace core
#endif
