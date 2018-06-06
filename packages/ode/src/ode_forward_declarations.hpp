
#ifndef ODE_FORWARD_DECLARATIONS_HPP_
#define ODE_FORWARD_DECLARATIONS_HPP_

#include "ode_ConfigDefs.hpp"


namespace ode {
  
  // template<typename state_type,
  // 	   typename rhs_type,
  // 	   typename scalar_type,
  // 	   typename state_resizer_fnctor_type,
  // 	   typename enable = void
  // 	   >
  // class rungeKutta4Stepper;


  template<typename state_type,
	   typename rhs_type,
	   typename scalar_type,
	   typename state_resizer_fnctor_type,
	   typename residual_policy_type,
	   typename enable = void
	   >
  class explicitEulerStepper;
  
  template<typename state_type,
	   typename rhs_type,
	   typename jacobian_type,
	   typename scalar_type,
	   typename state_resizer_fnctor_type,
	   typename residual_policy_type,
	   typename jacobian_policy_type,
	   typename nonlinearsolver_policy_type,
	   typename enable = void
	   >
  class implicitEulerStepper;

  
} // end namespace 

#endif
