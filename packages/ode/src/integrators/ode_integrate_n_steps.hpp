
#ifndef ODE_INTEGRATE_N_STEPS_HPP_
#define ODE_INTEGRATE_N_STEPS_HPP_

#include "ode_ConfigDefs.hpp"
#include "ode_integrate_n_steps_impl.hpp"


namespace ode{
  
template<typename stepper_type,
	 typename state_type,
	 typename collector_functor_type,
	 typename time_type,
	 typename std::enable_if<
	   !std::is_void<stepper_type>::value &&
	   !std::is_void<collector_functor_type>::value && 
	   !std::is_same<collector_functor_type,
			 voidCollector>::value
	   >::type * = nullptr
	 >
void integrateNSteps(stepper_type & stepper,
		     state_type & stateIn,
		     collector_functor_type & collector, 
		     time_type start_time,
		     time_type dt,
		     size_t num_steps)
{
  impl::integrateNStepsImpl<stepper_type, state_type,
			time_type, collector_functor_type>
    (stepper,stateIn,start_time, dt,num_steps, collector);
}
//----------------------------------------------------------------


template<typename stepper_type,
	 typename state_type,
	 typename time_type,
	 typename std::enable_if<
	   !std::is_void<stepper_type>::value
	   >::type * = nullptr
	 >
void integrateNSteps(stepper_type & stepper,
		     state_type & stateIn,
		     time_type start_time,
		     time_type dt,
		     size_t num_steps)
{
  impl::integrateNStepsImpl<stepper_type,
			    state_type,
			    time_type>
    (stepper, stateIn, start_time, dt, num_steps);
}

}//end namespace
#endif 
