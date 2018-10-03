
#ifndef ODE_FORWARD_DECLARATIONS_HPP_
#define ODE_FORWARD_DECLARATIONS_HPP_

#include "ode_ConfigDefs.hpp"

namespace rompp{
namespace ode{
  
// template<typename ode_state_type,
// 	 typename ode_residual_type,
// 	 typename model_type,
// 	 typename residual_policy_type = void,
// 	 typename enable = void
// 	 >
// class ExplicitEulerStepper;

template<ExplicitSteppersEnum whichone,
	 typename ode_state_type,
	 typename ... Args>
class ExplicitStepper;
  
template<typename ode_state_type,
	 typename ode_residual_type,
	 typename model_type,
	 typename residual_policy_type = void,
	 typename enable = void
	 >
class ExplicitRungeKutta4Stepper;
  
template<typename ode_state_type,
         typename ode_residual_type,
         typename ode_jacobian_type,
         typename model_type,
         typename residual_policy_type=void,
         typename jacobian_policy_type=void,
	 typename enable = void
         >
class ImplicitEulerStepper;
    
} // end namespace ode


/////////////////////////////////////
/////////////////////////////////////


namespace ode{
namespace impl {

template<typename ode_state_type,
	 typename ode_residual_type,
	 typename scalar_type,
	 typename model_type,	
	 typename residual_policy_type,
	 typename enable = void
	 >
class ExplicitEulerStepperImpl;

template<typename ode_state_type,
	 typename ode_residual_type,
	 typename scalar_type,
	 typename model_type,	
	 typename residual_policy_type,
	 typename enable = void
	 >
class ExplicitRungeKutta4StepperImpl;
  
template<typename ode_state_type,
	 typename ode_residual_type,
	 typename ode_jacobian_type,
	 typename scalar_type,
	 typename model_type,
	 typename residual_policy_type,
	 typename jacobian_policy_type,
	 typename enable = void
	 >
class ImplicitEulerStepperImpl;
  
}//end namespace impl
} // end namespace ode
}//end namespace rompp
#endif







// template<typename state_type,
//          typename residual_type,
//          typename jacobian_type,
//          typename scalar_type,
//          typename model_type,
// 	 typename time_type,
// 	 typename solver_policy_type,
// 	 typename aux_start_stepper_type,
//          typename residual_policy_type=void,
//          typename jacobian_policy_type=void,
// 	 typename enable = void
//          >
// class implicitBDF2Stepper;

// template<typename state_type,
//          typename residual_type,
//          typename jacobian_type,
//          typename scalar_type,
//          typename model_type,
// 	 typename time_type,
// 	 typename solver_policy_type,
// 	 typename aux_start_stepper_type,
//          typename residual_policy_type=void,
//          typename jacobian_policy_type=void,
// 	 typename enable = void
//          >
// class implicitBDF3Stepper;

// template<typename state_type,
//          typename residual_type,
//          typename jacobian_type,
//          typename scalar_type,
//          typename model_type,
// 	 typename time_type,
// 	 typename solver_policy_type,
//          typename residual_policy_type=void,
//          typename jacobian_policy_type=void,
// 	 typename enable = void
//          >
// class implicitAdamsMoulton1Stepper;


// template<typename state_type,
// 	 typename residual_type,
// 	 typename scalar_type,
// 	 typename model_type,	
// 	 typename residual_policy_type,
// 	 typename butcher_table_type,
// 	 typename enable = void
// 	 >
// class explicitAnyRungeKuttaStepperImpl;


// template<typename state_type,
//          typename residual_type,
//          typename jacobian_type,
//          typename scalar_type,
//          typename model_type,
// 	 typename solver_policy_type,
// 	 typename aux_start_stepper_type,
//          typename residual_policy_type,
//          typename jacobian_policy_type,
// 	 typename enable = void
//          >
// class implicitBDF2StepperImpl;

// template<typename state_type,
//          typename residual_type,
//          typename jacobian_type,
//          typename scalar_type,
//          typename model_type,
// 	 typename solver_policy_type,
// 	 typename aux_start_stepper_type,
//          typename residual_policy_type,
//          typename jacobian_policy_type,
// 	 typename enable = void
//          >
// class implicitBDF3StepperImpl;

// template<typename state_type,
//          typename residual_type,
//          typename jacobian_type,
//          typename scalar_type,
//          typename model_type,
// 	 typename solver_policy_type,
//          typename residual_policy_type,
//          typename jacobian_policy_type,
// 	 typename enable = void
//          >
// class implicitAdamsMoulton1StepperImpl;
