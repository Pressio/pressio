
#ifndef ODE_IMPLICIT_EULER_POLICIES_META_HPP_
#define ODE_IMPLICIT_EULER_POLICIES_META_HPP_

#include "ode_implicit_policies_meta.hpp"
#include "../standard/ode_implicit_euler_residual_standard_policy.hpp"
#include "../standard/ode_implicit_euler_jacobian_standard_policy.hpp"

namespace ode{
namespace meta {

//-----------------------------------------------------------------
// METAF FOR ADMISSIBLE IMPLICIT EULER RESIDUAL
//-----------------------------------------------------------------
template<typename policy_t, typename enable = void>
struct isLegitimateImplicitEulerResidualPolicy : std::false_type{};
  
template <typename policy_t>
struct isLegitimateImplicitEulerResidualPolicy<
  policy_t,
  typename std::enable_if<
    core::meta::publiclyInheritsFrom<
      policy_t,
      ode::policy::implicitResidualPolicyBase<policy_t, 1, 0>
      >::value 
    >::type
  > : std::true_type{};
  
  
//-----------------------------------------------------------------
// METAF FOR ADMISSIBLE IMPLICIT EULER JACOBIAN
//-----------------------------------------------------------------
template<typename policy_t, typename enable = void>
struct isLegitimateImplicitEulerJacobianPolicy
  : isLegitimateImplicitJacobianPolicy<policy_t>{};


//-----------------------------------------------------------------
// METAF TO CHECK RESIDUAL POLICY IS STANDARD
//-----------------------------------------------------------------
template<typename policy_t, typename enable = void>
struct isImplicitEulerResidualStandardPolicy
  : std::false_type{};

  
template <template <typename...> class policy_t,
	  typename state_type,
	  typename residual_type,
	  typename model_type,
	  typename time_type,
	  typename sizer_type>
struct isImplicitEulerResidualStandardPolicy<
  policy_t<state_type, residual_type, model_type, 
	   time_type, sizer_type>,
  typename std::enable_if<
    std::is_same<policy_t<state_type, residual_type, 
			  model_type, time_type, sizer_type>,
		 ode::policy::implicitEulerResidualStandardPolicy<
		   state_type, residual_type, 
		   model_type, time_type, sizer_type>
		 >::value
			  >::type
  > : std::true_type{};

  
//-----------------------------------------------------------------
// METAF TO CHECK JACOBIAN POLICY IS STANDARD
//-----------------------------------------------------------------
template<typename policy_t, typename enable = void>
struct isImplicitEulerJacobianStandardPolicy
  : std::false_type{};

  
template <template <typename...> class policy_t,
	  typename state_type,
	  typename jacobian_type,
	  typename model_type,
	  typename time_type,
	  typename sizer_type>
struct isImplicitEulerJacobianStandardPolicy<
  policy_t<state_type, jacobian_type, model_type,
	   time_type, sizer_type>,
  typename std::enable_if<
    std::is_same<policy_t<state_type, jacobian_type, 
			  model_type, time_type, sizer_type>,
		 ode::policy::implicitEulerJacobianStandardPolicy<
		   state_type, jacobian_type, 
		   model_type, time_type, sizer_type>
		 >::value
			  >::type
  > : std::true_type{};
//----------------------------------------------------------------

  
} // namespace meta
} // namespace core
#endif
