
#ifndef TIMEINTEGRATOR_STEPPER_EULER_HPP
#define TIMEINTEGRATOR_STEPPER_EULER_HPP

#include "timeIntegrator_explicit_stepper_base.hpp"

template<
  class state_type,
  class deriv_type = state_type,
  class scalar_type = typename timeIntegrator::details::defaultTypes::scalar_t,
  class time_type = typename timeIntegrator::details::defaultTypes::scalar_t
>
class euler
: public explicit_stepper_base<
  euler< state_type, deriv_type, scalar_type, time_type>,
  1, state_type, deriv_type, scalar_type, time_type>
{
public :
  using stepper_base_t = explicit_stepper_base<
  euler< state_type, deriv_type, scalar_type, time_type>,
  1, state_type, deriv_type, scalar_type, time_type>;

  // (de)constructors
  euler() : stepper_base_t(){}
  virtual ~euler(){}

  // methods
  template< class system_type>
  void step_impl( System /* system */,
		  const state_type &in,
		  const deriv_type & rhs,
		  state_type & out,
		  time_type /*t*/,
		  time_type dt )
  {
    cassert( in.size()==out.size() && in.size()==rhs.size() );
    
    // TODO: change this so that it uses the methods of the target data types
    // out = in + dt * rhs
    for (int i=0; i<in.size(); i++){
      out[i] = in[i] + dt*rhs[i];
    }
  }


  
}; //end class

  
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
