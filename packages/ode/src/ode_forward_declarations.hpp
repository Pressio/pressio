
#ifndef ODE_FORWARD_DECLARATIONS_HPP_
#define ODE_FORWARD_DECLARATIONS_HPP_

#include "ode_ConfigDefs.hpp"
//#include "core_forward_declarations.hpp"

namespace svd {

  // namespace details {
  //   template<typename T, typename enable = void>
  //   struct traits; // : core::details::traits<T,enable>{};
  //   template<typename T>
  //   struct traits<const T> : traits<T>{}; //core::details::traits<const T> {};
  // } // end namespace details
  
  // forward declaration 

  template<typename stepper_type,
	   typename functor_type, 
	   typename state_type,
	   typename observer_type,
	   typename enable = void
	   >
  void integrate_n_steps;


  template<typename state_type,
	   typename residual_type,
	   typename scalar_type,
	   typename enable = void
	   >
  class eulerStepper;

  
  
  
} // end namespace 

#endif
