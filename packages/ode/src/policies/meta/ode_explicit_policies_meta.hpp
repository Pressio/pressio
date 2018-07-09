
#ifndef ODE_EXPLICIT_POLICIES_META_HPP_
#define ODE_EXPLICIT_POLICIES_META_HPP_

#include "ode_policies_meta.hpp"
#include "../standard/ode_residual_standard_policy.hpp"


namespace ode{
namespace meta {
  
template<typename policy_t,
	 typename enable = void>
struct isLegitimateExplicitResidualPolicy 
: isLegitimateResidualPolicy<policy_t>{};

//-----------------------------------------------------------------
  
template<typename policy_t,
	 typename enable = void>
struct isExplicitEulerResidualStandardPolicy 
: isResidualStandardPolicy<policy_t>{};

//----------------------------------------------------------------

template<typename policy_t,
	 typename enable = void>
struct isExplicitRungeKutta4ResidualStandardPolicy
: isResidualStandardPolicy<policy_t>{};

//----------------------------------------------------------------
  

} // namespace meta
} // namespace core
#endif
