
#ifndef ODE_EULER_STEPPER_HPP_
#define ODE_EULER_STEPPER_HPP_

#include "ode_explicit_stepper_base.hpp"


template<typename state_type,
	 typename rhs_type,
	 typename scalar_type
	 >
class eulerStepper<state_type, rhs_type, scalar_type,
		   std::enable_if< !std::is_same<state_type,void>::value
				   >::type>
  : public explicit_stepper_base<eulerStepper<state_type,rhs_type,scalar_type> >
				 /*1,state_type,rhs_type,scalar_type*/
{
public :
  using stepper_type = eulerStepper<state_type,deriv_type,scalar_type>;
  using stepper_base_t = explicit_stepper_base<stepper_type>;
				/*,1,state_type,rhs_type,
				  scalar_type,ode::details::time_type*/

  // (de)constructors
  eulerStepper() : stepper_base_t(){}
  ~eulerStepper(){}

  // methods
  template< typename functor_type>
  void do_step_impl(functor_type & func,
		    const state_type &in,
		    state_type & out,
		    ode::details::time_type /*t*/,
		    ode::details::time_type dt )
  {
    cassert( in.size()==out.size() && in.size()==rhs.size() );

    functor(inout, R_, t);
    
    // TODO: if possible, change to use native operations of the target data types
    // out = in + dt * rhs
    for (int i=0; i<in.size(); i++){
      out[i] = in[i] + dt*rhs[i];
    }
  }

private:
    rhs_type R_;
  
}; //end class

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
