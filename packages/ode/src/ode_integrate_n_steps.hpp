
#ifndef ODE_INTEGRATE_N_STEPS_HPP_
#define ODE_INTEGRATE_N_STEPS_HPP_

#include "ode_ConfigDefs.hpp"
#include "ode_forward_declarations.hpp"
#include "ode_void_collector.hpp"

namespace ode{

  
  template<typename stepper_type,
  	   typename rhs_functor_type,	   
  	   typename state_type,
  	   typename time_type,
  	   typename n_steps_type,
  	   typename collector_type
  	   >
  typename std::enable_if<!std::is_void<stepper_type>::value &&
  			  !std::is_void<rhs_functor_type>::value &&
			  std::is_integral<n_steps_type>::value
  			  >::type
  integrateNStepsImpl(stepper_type & stepper,
  		      rhs_functor_type & rhs_functor,
  		      state_type & stateIn,
  		      time_type start_time,
  		      time_type dt,
  		      n_steps_type num_of_steps,
  		      typename std::conditional< std::is_same<collector_type, ode::voidCollector>::value,
		      ode::voidCollector, collector_type &>::type collector)
  {
    ode::details::time_type time = start_time;
    size_t step = 0;
    for( ; step < num_of_steps ; ++step)
    {
      // call collector/observer
      collector(step, time, stateIn);

      // do one step
      stepper.doStep( rhs_functor, stateIn, time, dt );

      // advance time
      time += + dt;
    }
    collector(step, time, stateIn);
  }




  /*+++++++++++++++++++++++++++++++++++++++++++++++++++++
  
    TODO: separate interface and impl

  +++++++++++++++++++++++++++++++++++++++++++++++++++++*/  

  template<typename stepper_type,
	   typename rhs_functor_type, 
	   typename state_type,
	   typename collector_functor_type
	   >
  typename std::enable_if<!std::is_same<stepper_type,void>::value &&
			  !std::is_same<rhs_functor_type,void>::value 
			  >::type
  integrateNSteps(stepper_type & stepper,
		  rhs_functor_type & rhs_functor,
		  state_type & stateIn,
		  collector_functor_type & collector,
		  ode::details::time_type start_time,
		  ode::details::time_type dt,
		  size_t num_steps)
  {
    integrateNStepsImpl<stepper_type, rhs_functor_type, state_type,
      			ode::details::time_type, size_t,
      			collector_functor_type>(stepper, rhs_functor, stateIn,
    			start_time, dt, num_steps, collector);
  }


  template<typename stepper_type,
	   typename rhs_functor_type,
	   typename state_type
	   >
  typename std::enable_if<!std::is_same<stepper_type,void>::value &&
			  !std::is_same<rhs_functor_type,void>::value 
			  >::type
  integrateNSteps(stepper_type & stepper,
		  rhs_functor_type & rhs_functor,
		  state_type & stateIn,
		  ode::details::time_type start_time,
		  ode::details::time_type dt,
		  size_t num_steps)
  {
    integrateNStepsImpl<stepper_type, rhs_functor_type, state_type,
      			ode::details::time_type, size_t,
			ode::voidCollector>(stepper, rhs_functor,stateIn,
					    start_time, dt, num_steps,
					    ode::voidCollector());
  }

  
}//end namespace

#endif 
