
#ifndef ODE_FORWARD_DECLARATIONS_HPP_
#define ODE_FORWARD_DECLARATIONS_HPP_

#include "ode_ConfigDefs.hpp"


namespace ode {

  namespace details {
    template<typename T, typename enable = void>
    struct traits; // : core::details::traits<T,enable>{};

    template<typename T>
    struct traits<const T> : traits<T>{}; //core::details::traits<const T> {};
  } // end namespace details
  

  // forward declaration 
  template<typename state_type,
	   typename rhs_type,
	   typename scalar_type,
	   typename state_resizer_fnctor_type, 
	   typename enable = void
	   >
  class eulerStepper;

  template<typename state_type,
	   typename rhs_type,
	   typename scalar_type,
	   typename state_resizer_fnctor_type,
	   typename enable = void
	   >
  class rungeKutta4Stepper;

  template<typename state_type,
	   typename rhs_type,
	   typename jacobian_type,
	   typename functor_type,
	   typename scalar_type,
	   typename state_resizer_fnctor_type,
	   typename enable = void
	   >
  class implicitEulerStepper;

  
} // end namespace 

#endif
