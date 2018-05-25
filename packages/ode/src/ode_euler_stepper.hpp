
#ifndef ODE_EULER_STEPPER_HPP_
#define ODE_EULER_STEPPER_HPP_

#include "ode_explicit_stepper_base.hpp"


namespace ode{

template<typename state_type,
	 typename rhs_type,
	 typename scalar_type,
	 typename state_resizer_fnctor_type
	 >
class eulerStepper<state_type, rhs_type, scalar_type, state_resizer_fnctor_type,
		   typename 
		   std::enable_if< !std::is_same<state_type,void>::value &&
				   core::meta::is_default_constructible<state_resizer_fnctor_type>::value
				   >::type
		   >
  : public explicit_stepper_base<eulerStepper<state_type,rhs_type,
					      scalar_type, state_resizer_fnctor_type> >
{
public :
  using stepper_t = eulerStepper<state_type,rhs_type,scalar_type,state_resizer_fnctor_type>;
  using stepper_base_t = explicit_stepper_base<stepper_t>;

  // (de)constructors
  eulerStepper() : stepper_base_t(){}
  ~eulerStepper(){}


  // methods
  template< typename functor_type>
  void doStepImpl(functor_type & functor,
		  state_type & y_inout,
		  ode::details::time_type t,
		  ode::details::time_type dt )
  {
    this->myResizer_(y_inout, RHS_);

    //eval RHS
    functor(y_inout, RHS_, t);
    
    // TODO: if possible, change to use native operations of the target data types
    // out = in + dt * rhs    
    for (decltype(y_inout.size()) i=0; i < y_inout.size(); i++){
      y_inout[i] += dt*RHS_[i];
    }
  }

private:
    rhs_type RHS_;
  
}; //end class


}//end namespace

#endif 




  
  // #ifndef DOXYGEN_SKIP
    // #else
    // typedef explicit_stepper_base< euler< ... > , ... > stepper_base_type;
    // #endif
    // typedef typename stepper_base_type::state_type state_type;
    // typedef typename stepper_base_type::value_type value_type;
    // typedef typename stepper_base_type::deriv_type deriv_type;
    // typedef typename stepper_base_type::time_type time_type;
    // typedef typename stepper_base_type::algebra_type algebra_type;
    // typedef typename stepper_base_type::operations_type operations_type;
    // typedef typename stepper_base_type::resizer_type resizer_type;

    // #ifndef DOXYGEN_SKIP
    // typedef typename stepper_base_type::stepper_type stepper_type;
    // typedef typename stepper_base_type::wrapped_state_type wrapped_state_type;
    // typedef typename stepper_base_type::wrapped_deriv_type wrapped_deriv_type;
    // #endif 
