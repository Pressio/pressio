
#ifndef ODE_INTEGRATORS_INTEGRATE_N_STEPS_IMPL_HPP_
#define ODE_INTEGRATORS_INTEGRATE_N_STEPS_IMPL_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../ode_forward_declarations.hpp"

namespace rompp{ namespace ode{ namespace impl{

struct voidCollector{
  template <typename step_t, typename scalar_t, typename state_t>
  void operator()(step_t step, scalar_t t, const state_t & y){
    //no op
  }
};


template<typename stepper_type,	typename state_type, typename time_type,
	 typename integral_type, bool is_implicit, typename solver_type>
void doStep(stepper_type & stepper, state_type & yIn,
	    time_type time, time_type dt, integral_type step,
	    core::meta::enable_if_t<is_implicit==true, solver_type> * solver){
  stepper(yIn, time, dt, step, *solver);
}

template<typename stepper_type,	typename state_type, typename time_type,
	 typename integral_type, bool is_implicit, typename solver_type = void>
void doStep(stepper_type & stepper, state_type & yIn,
	    time_type time, time_type dt, integral_type step,
	    core::meta::enable_if_t<is_implicit==false, solver_type> * solver = nullptr){
  stepper(yIn, time, dt, step);
}


template<typename stepper_type,	typename state_type, typename time_type,
	 typename integral_type, bool is_implicit,
	 typename collector_type = voidCollector,
	 typename solver_type = core::impl::empty>
void integrateNSteps(stepper_type   & stepper,
		     state_type     & yIn,
		     time_type	      start_time,
		     time_type	      dt,
		     integral_type    num_steps,
		     collector_type * collector = nullptr,
		     solver_type    * solver	= nullptr){

  time_type time = start_time;
  // call collector/observer at starting time
  if (collector) (*collector)(0, time, yIn);

#ifdef HAVE_TEUCHOS_TIMERS
  auto timer = Teuchos::TimeMonitor::getStackedTimer();
  timer->start("time loop");
#endif

  // time loop
  integral_type step = 1;
  for( ; step <= num_steps ; ++step){

#ifdef HAVE_TEUCHOS_TIMERS
    timer->start("time step");

    doStep<stepper_type, state_type,
	   time_type, integral_type,
	   is_implicit, solver_type>(stepper, yIn, time,
				     dt, step, solver);

    timer->stop("time step");
#endif

    // advance time
    time = start_time + static_cast<time_type>(step) * dt;

    // call collector/observer
    if (collector) (*collector)(step, time, yIn);
  }
#ifdef HAVE_TEUCHOS_TIMERS
  timer->stop("time loop");
#endif

}//edn integrateNSteps


}}}//end namespace rompp::ode::impl
#endif
