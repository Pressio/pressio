
#ifndef ODE_FORWARD_DECLARATIONS_HPP_
#define ODE_FORWARD_DECLARATIONS_HPP_

#include "ode_ConfigDefs.hpp"

namespace ode {

template<typename state_type,
	 typename residual_type,
	 typename scalar_type,
	 typename model_type,
	 typename time_type,
	 typename sizer_type,
	 typename residual_policy_type = void,
	 typename enable = void
	 >
class explicitEulerStepper;

template<typename state_type,
	 typename residual_type,
	 typename scalar_type,
	 typename model_type,
	 typename time_type,
	 typename sizer_type,
	 typename residual_policy_type = void,
	 typename enable = void
	 >
class explicitRungeKutta4Stepper;
  

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
// class implicitEulerStepper;


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
//          typename residual_policy_type=void,
//          typename jacobian_policy_type=void,
// 	 typename enable = void
//          >
// class implicitAdamsMoulton1Stepper;
  
  
} // end namespace ode


/////////////////////////////////////
/////////////////////////////////////


namespace ode {
namespace impl {

template<typename state_type,
	 typename residual_type,
	 typename scalar_type,
	 typename model_type,	
	 typename time_type,
	 typename sizer_type,
	 typename residual_policy_type,
	 typename enable = void
	 >
class explicitEulerStepperImpl;

template<typename state_type,
	 typename residual_type,
	 typename scalar_type,
	 typename model_type,	
	 typename time_type,
	 typename sizer_type,
	 typename residual_policy_type,
	 typename enable = void
	 >
class explicitRungeKutta4StepperImpl;
  

// template<typename state_type,
//          typename residual_type,
//          typename jacobian_type,
//          typename scalar_type,
//          typename model_type,
// 	 typename time_type,
// 	 typename solver_policy_type,
//          typename residual_policy_type,
//          typename jacobian_policy_type,
// 	 typename enable = void
//          >
// class implicitEulerStepperImpl;


// template<typename state_type,
//          typename residual_type,
//          typename jacobian_type,
//          typename scalar_type,
//          typename model_type,
// 	 typename time_type,
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
// 	 typename time_type,
// 	 typename solver_policy_type,
//          typename residual_policy_type,
//          typename jacobian_policy_type,
// 	 typename enable = void
//          >
// class implicitAdamsMoulton1StepperImpl;

  
}//end namespace impl  

} // end namespace ode
#endif
