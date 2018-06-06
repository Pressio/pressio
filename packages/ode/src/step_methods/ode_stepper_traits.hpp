

#ifndef ODE_STEPPER_TRAITS_HPP_
#define ODE_STEPPER_TRAITS_HPP_

#include "ode_forward_declarations.hpp"


namespace ode{
namespace details{
  
  template<typename state_type,
	   typename rhs_type,
	   typename scalar_type,
	   typename state_resizer_fnctor_type,
	   typename residual_policy_type
	   >
  struct traits< explicitEulerStepper<state_type, rhs_type, scalar_type,
				      state_resizer_fnctor_type, residual_policy_type> >
  {
    using order_t = unsigned int;
    using stepper_t = explicitEulerStepper<state_type,rhs_type,scalar_type,
				   state_resizer_fnctor_type, residual_policy_type>;
    using state_t =  state_type;
    using rhs_t = rhs_type;
    using scalar_t = scalar_type;
    using resizer_t = state_resizer_fnctor_type;
    using residual_policy_t = residual_policy_type;

    static constexpr order_t order_value = 1;    
  };

  
  // template<typename state_type,
  // 	   typename rhs_type,
  // 	   typename scalar_type,
  // 	   typename state_resizer_fnctor_type
  // 	   >
  // struct traits< rungeKutta4Stepper<state_type,rhs_type,scalar_type,state_resizer_fnctor_type> >
  // {
  //   using order_t = unsigned int;
  //   static constexpr order_t order_value = 4;

  //   using stepper_t = rungeKutta4Stepper<state_type,rhs_type,scalar_type,state_resizer_fnctor_type>;
  //   using state_t =  state_type;
  //   using rhs_t = rhs_type;
  //   using scalar_t = scalar_type;    
  //   using resizer_t = state_resizer_fnctor_type; 
  // };
  

  template<typename state_type,
	   typename rhs_type,
	   typename jacobian_type,
	   typename scalar_type,
	   typename state_resizer_fnctor_type,
	   typename residual_policy_type,
	   typename jacobian_policy_type,
	   typename nonlinearsolver_policy_type
	   >
  struct traits< implicitEulerStepper<state_type, rhs_type, jacobian_type, scalar_type,
				      state_resizer_fnctor_type, residual_policy_type,
				      jacobian_policy_type, nonlinearsolver_policy_type> >
  {
    using order_t = unsigned int;
    using stepper_t = implicitEulerStepper<state_type, rhs_type, jacobian_type,
					   scalar_type, state_resizer_fnctor_type,
					   residual_policy_type, jacobian_policy_type,
					   nonlinearsolver_policy_type>;
    using state_t =  state_type;
    using rhs_t = rhs_type;
    using jacobian_t =  jacobian_type;
    using scalar_t = scalar_type;    
    using resizer_t = state_resizer_fnctor_type;
    using residual_policy_t = residual_policy_type;
    using jacobian_policy_t = jacobian_policy_type;
    using nonlinearsolver_policy_t = nonlinearsolver_policy_type;

    static constexpr order_t order_value = 1;
  };


  
}//end namespace details
}//end namespace ode

#endif
