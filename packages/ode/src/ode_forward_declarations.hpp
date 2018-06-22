
#ifndef ODE_FORWARD_DECLARATIONS_HPP_
#define ODE_FORWARD_DECLARATIONS_HPP_

#include "ode_ConfigDefs.hpp"
// standard policies for explicit methods
#include "./policies/ode_explicit_euler_standard_policy.hpp"
#include "./policies/ode_explicit_runge_kutta4_standard_policy.hpp"
// standard policies for implicit methods
#include "./policies/ode_implicit_euler_residual_standard_policy.hpp"
#include "./policies/ode_implicit_euler_jacobian_standard_policy.hpp"


namespace ode {

    
template<typename state_type,
	 typename residual_type,
	 typename scalar_type,
	 typename state_resizer_fnctor_type,
	 typename model_type,
	 typename time_type,
	 typename residual_policy_type =
	   ode::policy::explicitEulerStandardResidual<
	   state_type,residual_type,model_type,time_type>,
	 typename enable = void
	 >
class explicitEulerStepper;

  

template<typename state_type,
	 typename residual_type,
	 typename scalar_type,
	 typename state_resizer_fnctor_type,
	 typename model_type,
	 typename time_type,
	 typename residual_policy_type =
	   ode::policy::explicitRungeKutta4StandardResidual<
	   state_type,residual_type,model_type,time_type>,
	 typename enable = void
	 >
class explicitRungeKutta4Stepper;

  
template<typename state_type,
         typename residual_type,
         typename jacobian_type,
         typename scalar_type,
         typename state_resizer_fnctor_type,
         typename model_type,
	 typename time_type,
	 typename solver_policy_type,
         typename residual_policy_type = 
           ode::policy::implicitEulerStandardResidual<
	   state_type,residual_type,model_type,time_type>,
         typename jacobian_policy_type =
	   ode::policy::implicitEulerStandardJacobian<
	   state_type,jacobian_type,model_type,time_type>,
	 typename enable = void
         >
class implicitEulerStepper;


  
} // end namespace 

#endif
