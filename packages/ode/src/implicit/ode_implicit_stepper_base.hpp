
#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_BASE_IMPLICIT_STEPPER_BASE_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_BASE_IMPLICIT_STEPPER_BASE_HPP_

namespace pressio{ namespace ode{ namespace implicitmethods{

template<typename derived_type>
class StepperBase
{
  // friend the derived so that it can access private constructors
  friend derived_type;

public:
  types::stepper_order_t order() const{
    return static_cast<derived_type&>(*this).order();
  }

  template <typename state_t, typename scalar_t, typename solver_type>
  void doStep(state_t & odeState,
		  const scalar_t & time,
		  const scalar_t & dt,
		  const types::step_t & step,
		  solver_type & solver)
  {
    static_cast<derived_type &>(*this).template doStep<solver_type>(odeState, time, dt, step, solver);
  }

  template <typename state_t, typename scalar_t, 
  typename solver_type, typename guess_cb_t>
  void doStep(state_t & odeState,
		  const scalar_t & time,
		  const scalar_t & dt,
		  const types::step_t & step,
		  solver_type & solver,
		  guess_cb_t && guesser)
  {
    static_cast<derived_type &>(*this).template doStep<
      solver_type, guess_cb_t>(odeState, time, dt, step, solver, std::forward<guess_cb_t>(guesser));
  }

private:
  StepperBase() = default;
  ~StepperBase() = default;
  StepperBase(const StepperBase & other)  = delete;
  StepperBase & operator=(const StepperBase & other)  = delete;
  StepperBase(StepperBase && other)  = delete;
  StepperBase & operator=(StepperBase && other)  = delete;

};//end class

}}}//end namespace pressio::ode::implicitmethods
#endif
