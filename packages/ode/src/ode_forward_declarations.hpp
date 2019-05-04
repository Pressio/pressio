
#ifndef ODE_FORWARD_DECLARATIONS_HPP_
#define ODE_FORWARD_DECLARATIONS_HPP_

#include "ode_ConfigDefs.hpp"

namespace rompp{ namespace ode{

///begin namespace policy
namespace policy{

template<typename state_type,
	 typename model_type,
	 typename residual_type = state_type,
	 typename enable = void>
class ExplicitResidualStandardPolicy;

template<typename state_type,
	 typename model_type,
	 typename residual_type = state_type,
	 typename enable = void>
class ImplicitResidualStandardPolicy;

template<typename state_type,
	 typename model_type,
	 typename jacobian_type,
	 typename enable = void>
class ImplicitJacobianStandardPolicy;

}//end namespace policy
//-----------------------------------

template<ExplicitEnum whichone,
	 typename ode_state_type,
	 typename model_type,
	 typename ode_residual_type,
	 typename ...Args
	 >
class ExplicitStepper;

template<ImplicitEnum whichone,
	 typename ode_state_type,
	 typename ode_residual_type,
	 typename ode_jacobian_type,
	 typename model_type,
	 typename ...Args>
class ImplicitStepper;
//-----------------------------------

///begin namespace impl
namespace impl {

template<typename scalar_type,
	 typename ode_state_type,
	 typename model_type,
	 typename ode_residual_type,
	 typename residual_policy_type,
	 typename ops_t,
	 typename enable = void
	 >
class ExplicitEulerStepperImpl;

template<typename scalar_type,
	 typename ode_state_type,
	 typename model_type,
	 typename ode_residual_type,
	 typename residual_policy_type,
	 typename enable = void
	 >
class ExplicitRungeKutta4StepperImpl;

template<typename ode_state_type,
	 typename ode_residual_type,
	 typename ode_jacobian_type,
	 typename model_type,
	 typename residual_policy_type,
	 typename jacobian_policy_type,
	 typename enable = void
	 >
class ImplicitEulerStepperImpl;

template<typename ode_state_type,
	 typename ode_residual_type,
	 typename ode_jacobian_type,
	 typename model_type,
	 typename aux_stepper_type,
	 typename residual_policy_type,
	 typename jacobian_policy_type,
	 typename enable = void
	 >
class ImplicitBDF2StepperImpl;

}//end namespace impl

}} // end namespace rompp::ode
#endif
