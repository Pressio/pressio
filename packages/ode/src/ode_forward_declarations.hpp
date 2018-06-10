
#ifndef ODE_FORWARD_DECLARATIONS_HPP_
#define ODE_FORWARD_DECLARATIONS_HPP_

#include "ode_ConfigDefs.hpp"
#include "step_methods/policy/ode_explicit_euler_standard_policy.hpp"
#include "step_methods/policy/ode_implicit_euler_standard_policy.hpp"
#include "step_methods/policy/ode_explicit_runge_kutta4_standard_policy.hpp"


namespace ode {
  
  template<typename state_type,
	   typename residual_type,
	   typename scalar_type,
	   typename state_resizer_fnctor_type,
	   typename model_type,
	   typename residual_policy_type =
	   ode::policy::explicitRungeKutta4StandardResidual<state_type, residual_type,
							    model_type,
							    details::time_type>,
	   typename enable = void
	   >
  class explicitRungeKutta4Stepper;


  
  template<typename state_type,
	   typename residual_type,
	   typename scalar_type,
	   typename state_resizer_fnctor_type,
	   typename model_type,
	   typename residual_policy_type =
	   ode::policy::explicitEulerStandardResidual<state_type, residual_type,
						      model_type,
						      details::time_type>,
	   typename enable = void
	   >
  class explicitEulerStepper;

  
  
  template<typename state_type,
	   typename rhs_type,
	   typename jacobian_type,
	   typename scalar_type,
	   typename state_resizer_fnctor_type,
	   typename model_type,
	   typename residual_policy_type = 
	   ode::policy::implicitEulerStandardResidual<state_type,rhs_type,
						      model_type,details::time_type>,
	   typename jacobian_policy_type =
	            ode::policy::implicitEulerStandardJacobian<state_type,jacobian_type,
							       model_type,details::time_type>,
	   typename nonlinearsolver_policy_type =
	            void /*should be newton raphson or similar*/,
	   typename enable = void
	   >
  class implicitEulerStepper;

  
} // end namespace 

#endif
