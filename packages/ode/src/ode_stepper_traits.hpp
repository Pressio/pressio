

#ifndef ODE_STEPPER_TRAITS_HPP_
#define ODE_STEPPER_TRAITS_HPP_

#include "ode_forward_declarations.hpp"


namespace ode{
namespace details{
  
  template<typename state_type,
	   typename rhs_type,
	   typename scalar_type,
	   typename state_resizer_fnctor_type
	   >
  struct traits< eulerStepper<state_type, rhs_type, scalar_type, state_resizer_fnctor_type> >
  {
    using order_t = unsigned int;
    static constexpr order_t order_value = 1;

    using stepper_t = eulerStepper<state_type,rhs_type,scalar_type,state_resizer_fnctor_type>;
    using state_t =  state_type;
    using rhs_t = rhs_type;
    using scalar_t = scalar_type;    
    using resizer_t = state_resizer_fnctor_type; 
  };

  
  template<typename state_type,
	   typename rhs_type,
	   typename scalar_type,
	   typename state_resizer_fnctor_type
	   >
  struct traits< rungeKutta4Stepper<state_type,rhs_type,scalar_type,state_resizer_fnctor_type> >
  {
    using order_t = unsigned int;
    static constexpr order_t order_value = 4;

    using stepper_t = rungeKutta4Stepper<state_type,rhs_type,scalar_type,state_resizer_fnctor_type>;
    using state_t =  state_type;
    using rhs_t = rhs_type;
    using scalar_t = scalar_type;    
    using resizer_t = state_resizer_fnctor_type; 
  };
  

  template<typename state_type,
	   typename rhs_type,
	   typename jacobian_type,
	   typename functor_type,
	   typename scalar_type,
	   typename state_resizer_fnctor_type
	   >
  struct traits< implicitEulerStepper<state_type, rhs_type, jacobian_type, functor_type,
				      scalar_type, state_resizer_fnctor_type> >
  {
    using order_t = unsigned int;
    static constexpr order_t order_value = 1;

    using stepper_t = implicitEulerStepper<state_type, rhs_type, jacobian_type, functor_type,
					   scalar_type, state_resizer_fnctor_type>;
    using state_t =  state_type;
    using rhs_t = rhs_type;
    using jacobian_t =  jacobian_type;
    using scalar_t = scalar_type;    
    using resizer_t = state_resizer_fnctor_type;
    using functor_t = functor_type;
    //    using system_solver_t = system_solver_type;
  };


  
}//end namespace details
}//end namespace ode

#endif
