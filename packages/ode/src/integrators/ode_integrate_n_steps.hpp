
#ifndef ODE_INTEGRATE_N_STEPS_HPP_
#define ODE_INTEGRATE_N_STEPS_HPP_

#include "ode_ConfigDefs.hpp"
#include "../meta/ode_meta.hpp"
#include "./impl/ode_integrate_n_steps_impl.hpp"
#include "../steppers/implicit_steppers/ode_implicit_stepper_traits.hpp"
#include "../steppers/explicit_steppers/ode_explicit_stepper_traits.hpp"

namespace ode{
  
/* enable if: 
   (1) advancing the full state
   (2) collector is passed
*/
template<typename stepper_type, typename state_type,
	 typename time_type, typename integral_type,
	 typename collector_type,
	 typename std::enable_if<
	   ode::meta::isLegitimateCollector<collector_type, integral_type,
					    time_type, state_type>::value &&
	   ode::details::traits<typename
				stepper_type::base_t>::advanceIncrement==false
	   >::type * = nullptr >
void integrateNSteps(stepper_type & stepper,
		     state_type & yIn,
		     time_type start_time,
		     time_type dt,
		     integral_type num_steps,
		     collector_type & collector)
{
  impl::integrateNStepsImpl(stepper, yIn, start_time,
			    dt, num_steps, collector);
}  
//----------------------------------------------------------------

  
  
/* enable if: 
   (1) advancing the increment wrt y0
   (2) collector is passed
*/
template<typename stepper_type, typename state_type,
	 typename time_type, typename integral_type,
	 typename collector_type,
	 typename std::enable_if<
	   ode::meta::isLegitimateCollector<collector_type, integral_type,
					    time_type, state_type>::value &&
	   ode::details::traits<typename stepper_type::base_t>::advanceIncrement
	   >::type * = nullptr>
void integrateNSteps(stepper_type & stepper, state_type & yIn,
		     time_type start_time, time_type dt,
		     integral_type num_steps, collector_type & collector)
{
  using sizer_t = typename ode::details::traits<typename stepper_type::base_t>::sizer_t;
  auto stSz = sizer_t::getSize(yIn);
  state_type y(stSz);
  for (decltype(stSz) i=0; i<stSz; i++)
    y[i] = 0.0;

  impl::integrateNStepsImpl(stepper, y, start_time,
  			    dt, num_steps, collector);

  // yIn should contain the full final solution on exit
  for (decltype(stSz) i=0; i<stSz; i++)
    yIn[i] += y[i];
}
//----------------------------------------------------------------



  

// template<typename stepper_type,
// 	 typename state_type,
// 	 typename time_type,
// 	 typename integral_type,
// 	 typename std::enable_if<
// 	   !std::is_void<stepper_type>::value
// 	   >::type * = nullptr
// 	 >
// void integrateNSteps(stepper_type & stepper,
// 		     state_type & stateIn,
// 		     time_type start_time,
// 		     time_type dt,
// 		     integral_type num_steps)
// {
//   time_type time = start_time;
//   integral_type step = 0;
//   for( ; step < num_steps ; ++step)
//   {
//     // do one step
//     stepper.doStep(stateIn, time, dt);
//     // advance time: mulitply (vs adding) benefits roundoff
//     time = start_time + static_cast<time_type>(step) * dt;  
//   }
// }

}//end namespace
#endif 
