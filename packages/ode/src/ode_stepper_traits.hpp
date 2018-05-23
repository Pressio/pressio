

#ifndef ODE_STEPPER_TRAITS_HPP_
#define ODE_STEPPER_TRAITS_HPP_

#include "ode_forward_declarations.hpp"


namespace ode{
namespace details{

  template<typename state_type,
	   typename residual_type,
	   typename scalar_type
	   >
  struct traits< eulerStepper<state_type, residual_type, scalar_type> >
  {
    using order_t = unsigned int;
    static constexpr order_t order_value = 1;

    using stepper_t = eulerStepper<state_type, residual_type, scalar_type>;
    using state_t = typename state_type;
    using residual_t = typename residual_type;
    using scalar_t = typename scalar_type;
    
  };
  
  
}//End namespace details
}//end namespace ode

#endif
