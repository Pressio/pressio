
#ifndef ODE_INTEGRATORS_INTEGRATE_N_STEPS_IMPL_HPP_
#define ODE_INTEGRATORS_INTEGRATE_N_STEPS_IMPL_HPP_

#include "../ode_ConfigDefs.hpp"
#include "../ode_forward_declarations.hpp"

namespace rompp{ namespace ode{ namespace impl{

template< typename solver_type, typename guesser_cb_t>
struct DoStepMixin{
  template <typename time_type,
	    typename integral_type,
	    typename state_type,
	    typename stepper_type>
  void operator()(time_type time,
		  time_type dt,
		  integral_type step,
		  state_type & yIn,
		  stepper_type & stepper,
		  solver_type & solver,
		  guesser_cb_t && guessCb){
#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("time step");
#endif

    stepper(yIn, time, dt, step, solver,
	    std::forward<guesser_cb_t>(guessCb));

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("time step");
#endif
  }//end ()
};


template<>
struct DoStepMixin<core::impl::empty, core::impl::empty>{
  template <typename time_type,
	    typename integral_type,
	    typename state_type,
	    typename stepper_type>
  void operator()(time_type time,
		  time_type dt,
		  integral_type step,
		  state_type & yIn,
		  stepper_type & stepper){
#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("time step");
#endif

    stepper(yIn, time, dt, step);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("time step");
#endif
  }
};


template<typename solver_type>
struct DoStepMixin<solver_type, core::impl::empty>{
  template <typename time_type,
	    typename integral_type,
	    typename state_type,
	    typename stepper_type>
  void operator()(time_type time,
		  time_type dt,
		  integral_type step,
		  state_type & yIn,
		  stepper_type & stepper,
		  solver_type & solver){
#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("time step");
#endif

    stepper(yIn, time, dt, step, solver);

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("time step");
#endif
  }//end ()
};
//-------------------------------------------------------


/*
 * A valid collector object is passed by user
 * to take snapshots
 */
template <typename collector_type, typename DoStepMixin_t>
struct AdvancerMixin{

  template <typename integral_type,
	    typename time_type,
	    typename state_type,
	    typename ... Args>
  void operator()(integral_type    num_steps,
		  time_type	   start_time,
		  time_type	   dt,
		  state_type &	   yIn,
		  collector_type & collector,
		  Args && ...	   args)
  {
    DoStepMixin_t stepImpl;

#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("time loop");
#endif

    // time variable
    time_type time = start_time;
    // pass initial condition to collector object
    collector(0, time, yIn);

    integral_type step = 1;
    ::rompp::core::io::print_stdout("\nstarting time loop","\n");
    for( ; step <= num_steps ; ++step)
    {
      #ifdef DEBUG_PRINT
      auto fmt = core::io::bg_grey() + core::io::bold() + core::io::red();
      auto reset = core::io::reset();
      ::rompp::core::io::print_stdout(fmt, "doing time step =",
				      step, reset, "\n");
      #endif

      stepImpl(time, dt, step, yIn, std::forward<Args>(args)...);
      time = start_time + static_cast<time_type>(step) * dt;
      collector(step, time, yIn);
    }
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("time loop");
#endif
  }//end ()
};



/*
 * No collector object is passed by user
 */
template <typename DoStepMixin_t>
struct AdvancerMixin<core::impl::empty, DoStepMixin_t>{

  template <typename integral_type,
	    typename time_type,
	    typename ... Args>
  void operator()(integral_type num_steps,
		  time_type	start_time,
		  time_type	dt,
		  Args && ...	args)
  {
    DoStepMixin_t stepImpl;

#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("time loop");
#endif

    time_type time = start_time;
    integral_type step = 1;

    ::rompp::core::io::print_stdout("\nstarting time loop","\n");
    for( ; step <= num_steps ; ++step)
    {
      #ifdef DEBUG_PRINT
      auto fmt = core::io::bg_grey() + core::io::bold() + core::io::red();
      auto reset = core::io::reset();
      ::rompp::core::io::print_stdout(fmt, "doing time step =",
				      step, reset, "\n");
      #endif

      stepImpl(time, dt, step, std::forward<Args>(args)...);
      time = start_time + static_cast<time_type>(step) * dt;
    }

#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("time loop");
#endif

  }//end ()
};


}}}//end namespace rompp::ode::impl
#endif
