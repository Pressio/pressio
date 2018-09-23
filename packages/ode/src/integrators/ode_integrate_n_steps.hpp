
#ifndef ODE_INTEGRATORS_INTEGRATE_N_STEPS_HPP_
#define ODE_INTEGRATORS_INTEGRATE_N_STEPS_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../meta/ode_meta.hpp"

namespace rompp{
namespace ode{
  

//----------------------------------------------------------------
//----------------------------------------------------------------
//----------------------------------------------------------------
//                 FOR EXPLICIT METHODS
//----------------------------------------------------------------
//----------------------------------------------------------------
//----------------------------------------------------------------

 template<typename stepper_type,
	 typename state_type,
	 typename time_type,
	 typename integral_type,
	 typename collector_type,
	 typename std::enable_if<
	   ode::meta::isLegitimateCollector<collector_type,
					    integral_type,
					    time_type,
					    state_type>::value &&
	   std::is_integral<integral_type>::value &&
	   details::traits<typename
	     stepper_type::base_t>::is_explicit
	   >::type * = nullptr>
void integrateNSteps(stepper_type & stepper,
		     state_type & yIn,
		     time_type start_time,
		     time_type dt,
		     integral_type num_steps,
		     collector_type & collector)
{
  time_type time = start_time;
  integral_type step = 1;

  // time loop
  for( ; step <= num_steps ; ++step)
  {
    // call collector/observer 
    collector(step, time, yIn);
    // do one step
    stepper.doStep(yIn, time, dt, step);
    // advance time: mulitply (vs adding) benefits roundoff
    time = start_time + static_cast<time_type>(step) * dt;  
  }
  collector(step, time, yIn);
}  

//----------------------------------------------------------------

template<typename stepper_type,
	 typename state_type,
	 typename time_type,
	 typename integral_type,
	 typename std::enable_if<
	   details::traits<
	     typename stepper_type::base_t>::is_explicit
	   >::type * = nullptr>
void integrateNSteps(stepper_type & stepper,
		     state_type & yIn,
		     time_type start_time,
		     time_type dt,
		     integral_type num_steps)
{
  time_type time = start_time;
  integral_type step = 1;

  // time loop
  for( ; step <= num_steps ; ++step)
  {
    // do one step
    stepper.doStep(yIn, time, dt, step);
    // advance time: mulitply (vs adding) benefits roundoff
    time = start_time + static_cast<time_type>(step) * dt;  
  }
}  


  
//----------------------------------------------------------------
//----------------------------------------------------------------
//----------------------------------------------------------------
//                 FOR IMPLICIT METHODS
//----------------------------------------------------------------
//----------------------------------------------------------------
//----------------------------------------------------------------

template<typename stepper_type,
	 typename state_type,
	 typename time_type,
	 typename integral_type,
	 typename collector_type,
	 typename solver_type,
	 typename std::enable_if<
	   ode::meta::isLegitimateCollector<collector_type,
					    integral_type,
					    time_type,
					    state_type>::value &&
	   std::is_integral<integral_type>::value &&
	   details::traits<
	     typename stepper_type::base_t>::is_implicit
	   >::type * = nullptr>
void integrateNSteps(stepper_type & stepper,
		     state_type & yIn,
		     time_type start_time,
		     time_type dt,
		     integral_type num_steps,
		     collector_type & collector,
		     solver_type & solver)
{
  time_type time = start_time;
  integral_type step = 1;

  // time loop
  for( ; step <= num_steps ; ++step)
  {
    // call collector/observer 
    collector(step, time, yIn);
    // do one step
    stepper.doStep(yIn, time, dt, step, solver);
    // advance time: mulitply (vs adding) benefits roundoff
    time = start_time + static_cast<time_type>(step) * dt;  
  }
  collector(step, time, yIn);
}  

//----------------------------------------------------------------

template<typename stepper_type,
	 typename state_type,
	 typename time_type,
	 typename integral_type,
	 typename solver_type,
	 typename std::enable_if<
	   details::traits<
	     typename stepper_type::base_t>::is_implicit
	   >::type * = nullptr>
void integrateNSteps(stepper_type & stepper,
		     state_type & yIn,
		     time_type start_time,
		     time_type dt,
		     integral_type num_steps,
		     solver_type & solver)
{
  time_type time = start_time;
  integral_type step = 1;

  // time loop
  for( ; step <= num_steps ; ++step)
  {
    // do one step
    stepper.doStep(yIn, time, dt, step, solver);
    // advance time: mulitply (vs adding) benefits roundoff
    time = start_time + static_cast<time_type>(step) * dt;  
  }
}  

   
}//end namespace
}//end namespace rompp
#endif 
