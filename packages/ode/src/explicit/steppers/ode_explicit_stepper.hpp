
#ifndef ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_HPP_
#define ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_HPP_

#include "./impl/ode_explicit_euler_stepper_impl.hpp"
#include "ode_explicit_stepper_traits.hpp"

namespace rompp{ namespace ode{

// policy is STANDARD
template<ExplicitEnum whichone,
	 typename ode_state_type,
	 typename model_type,
	 typename ode_residual_type>
class ExplicitStepper<whichone, ode_state_type,
		      model_type, ode_residual_type, void>
  : public ExplicitStepperBase<
  ExplicitStepper<whichone, ode_state_type, model_type,
		  ode_residual_type, void> >{

  using this_t	 = ExplicitStepper<whichone, ode_state_type,
				   model_type, ode_residual_type, void>;
  using mytraits	 = details::traits<this_t>;
  using base_t	 = ExplicitStepperBase<this_t>;
  friend base_t;
  using scalar_type  = typename mytraits::scalar_t;
  using pol_t	 = typename mytraits::residual_policy_t;
  using impl_class_t = typename details::traits<this_t>::impl_t;

  impl_class_t myImpl_ = {};

public:
  ExplicitStepper() = delete;

  ~ExplicitStepper() = default;

  ExplicitStepper(const model_type & model,
		  ode_state_type const & y0,
		  ode_residual_type const & r0)
    : myImpl_(model, pol_t(), y0, r0){}


  template<typename ... Args>
  void operator()(Args && ... args){
    myImpl_.doStep( std::forward<Args>(args)... );
  }

  // template<typename step_t>
  // void operator()(ode_state_type & y, scalar_type t,
  // 		      scalar_type dt, step_t step){
  // 	myImpl_.doStep(y, t, dt, step);
  // }

};//end class


// policy is user-defined
template<ExplicitEnum whichone,
	 typename ode_state_type,
	 typename model_type,
	 typename ode_residual_type,
	 typename residual_policy_type>
class ExplicitStepper<whichone, ode_state_type,
		      model_type, ode_residual_type,
		      residual_policy_type>
  : public ExplicitStepperBase<
  ExplicitStepper<whichone, ode_state_type, model_type,
		  ode_residual_type, residual_policy_type> >{

  using this_t	 = ExplicitStepper<whichone, ode_state_type,
				   model_type, ode_residual_type,
				   residual_policy_type>;
  using mytraits	 = details::traits<this_t>;
  using base_t	 = ExplicitStepperBase<this_t>;
  friend base_t;
  using scalar_type  = typename mytraits::scalar_t;
  using impl_class_t = typename details::traits<this_t>::impl_t;

  impl_class_t myImpl_ = {};

public:
  ExplicitStepper(const model_type & model,
		  const residual_policy_type & policyObj,
		  ode_state_type const & y0,
		  ode_residual_type const & r0)
    : myImpl_(model, policyObj, y0, r0){}

  // when state and residual type are the same
  template <typename s_t = ode_state_type,
	    typename r_t = ode_residual_type,
	    core::meta::enable_if_t<
	      std::is_same<s_t, r_t>::value
	      > * = nullptr>
  ExplicitStepper(const model_type & model,
		  const residual_policy_type & policyObj,
		  s_t const & y0)
    : myImpl_(model, policyObj, y0, y0){}

  ExplicitStepper() = delete;
  ~ExplicitStepper() = default;

  template<typename step_t>
  void operator()(ode_state_type & y, scalar_type t,
		  scalar_type dt, step_t step){
    myImpl_.doStep(y, t, dt, step);
  }
};//end class

}} // end namespace rompp::ode
#endif
