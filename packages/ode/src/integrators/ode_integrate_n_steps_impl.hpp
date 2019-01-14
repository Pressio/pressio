
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


template <typename collector_type, typename DoStepMixin_t>
struct AdvancerMixin{

  template <typename integral_type,
	    typename time_type,
	    typename state_type,
	    typename ... Args>
  void operator()(integral_type num_steps,
		  time_type start_time,
		  time_type dt,
		  state_type & yIn,
		  collector_type & collector,
		  Args && ... args)
  {
    DoStepMixin_t stepImpl;

#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("time loop");
#endif
    time_type time = start_time;
    collector(0, time, yIn);

    integral_type step = 1;
    for( ; step <= num_steps ; ++step){
      stepImpl(time, dt, step, yIn, std::forward<Args>(args)...);
      time = start_time + static_cast<time_type>(step) * dt;
      collector(step, time, yIn);
    }
#ifdef HAVE_TEUCHOS_TIMERS
    timer->stop("time loop");
#endif
  }//end ()
};


template <typename DoStepMixin_t>
struct AdvancerMixin<core::impl::empty, DoStepMixin_t>{

  template <typename integral_type,
	    typename time_type,
	    typename ... Args>
  void operator()(integral_type num_steps,
		  time_type start_time,
		  time_type dt,
		  Args && ... args)
  {
    DoStepMixin_t stepImpl;

#ifdef HAVE_TEUCHOS_TIMERS
    auto timer = Teuchos::TimeMonitor::getStackedTimer();
    timer->start("time loop");
#endif

    time_type time = start_time;
    integral_type step = 1;

#ifdef DEBUG_PRINT
    ::rompp::core::debug::print("\nstarting time loop","\n");
#endif
    for( ; step <= num_steps ; ++step)
    {
#ifdef DEBUG_PRINT
      ::rompp::core::debug::print("\n--------------------","\n");
      ::rompp::core::debug::print("doing time step =", step, "\n");
      ::rompp::core::debug::print("---------------------","\n");
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











// #ifdef HAVE_TEUCHOS_TIMERS
//   auto timer = Teuchos::TimeMonitor::getStackedTimer();
//   timer->start("time loop");
// #endif

//   // time loop
//   integral_type step = 1;
//   for( ; step <= num_steps ; ++step){

// #ifdef HAVE_TEUCHOS_TIMERS
//     timer->start("time step");

//     doStep<stepper_type, state_type,
// 	   time_type, integral_type,
// 	   is_implicit, solver_type>(stepper, yIn, time,
// 				     dt, step, solver);

//     timer->stop("time step");
// #endif

//     // advance time
//     time = start_time + static_cast<time_type>(step) * dt;

//     // call collector/observer
//     if (collector) (*collector)(step, time, yIn);
//   }
// #ifdef HAVE_TEUCHOS_TIMERS
//   timer->stop("time loop");
// #endif
