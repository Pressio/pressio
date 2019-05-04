
#ifndef ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_HPP_
#define ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_HPP_

#include "./impl/ode_explicit_euler_stepper_impl.hpp"
#include "ode_explicit_stepper_traits.hpp"

namespace rompp{ namespace ode{

template<
  ExplicitEnum whichone,
  typename ode_state_type,
  typename model_type,
  typename ode_residual_type,
  typename ...Args
  >
class ExplicitStepper
  : public ExplicitStepperBase<
  ExplicitStepper<
    whichone, ode_state_type,
    model_type, ode_residual_type,
    Args...
    >
  >
{

  using this_t		= ExplicitStepper
    <whichone, ode_state_type, model_type, ode_residual_type, Args...>;
  using base_t		= ExplicitStepperBase<this_t>;
  friend base_t;

  using mytraits	= details::traits<this_t>;
  using scalar_type	= typename mytraits::scalar_t;
  using standard_res_policy_t = typename mytraits::standard_res_policy_t;
  using res_policy_t	= typename mytraits::residual_policy_t;
  using impl_class_t	= typename mytraits::impl_t;

  impl_class_t myImpl_ = {};

public:
  ExplicitStepper()  = delete;
  ~ExplicitStepper() = default;

  // only enable if the residual policy is standard
  template <typename T = ode_state_type,
	    ::rompp::mpl::enable_if_t<
	      mpl::is_same<
		standard_res_policy_t, res_policy_t
		>::value
	      > * = nullptr>
  ExplicitStepper(T const		  & y0,
		  const model_type	  & model,
		  ode_residual_type const & r0)
    : myImpl_(model, res_policy_t(), y0, r0){}


  // this is enabled all the time
  ExplicitStepper(ode_state_type const	  & y0,
		  const model_type	  & model,
		  const res_policy_t	  & policyObj,
		  ode_residual_type const & r0)
    : myImpl_(model, policyObj, y0, r0){}


  template<typename ... Args2>
  void operator()(Args2 && ... args){
    myImpl_.doStep( std::forward<Args2>(args)... );
  }

};//end class

}} // end namespace rompp::ode
#endif
