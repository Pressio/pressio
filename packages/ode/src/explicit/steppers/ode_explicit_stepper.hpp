
#ifndef ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_HPP_
#define ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_HPP_

#include "./impl/ode_explicit_euler_stepper_impl.hpp"
#include "ode_explicit_stepper_helper_info.hpp"

namespace rompp{ namespace ode{ 
      
template<ExplicitSteppersEnum whichone, typename... Args>
class ExplicitStepper
  : public impl::stepper_helper_info<whichone, Args...>::base_impl_type{

  using info_t = impl::stepper_helper_info<whichone, Args...>;
  using state_t = typename info_t::state_type;
  using model_t = typename info_t::model_type;
  using res_t = typename info_t::res_type;
  using pol_t = typename info_t::residual_policy_type;
  using pol_std_t = typename info_t::res_std_pol_type;
  using base_impl_t = typename info_t::base_impl_type;

//this needs to be public, it is detected by integrators
public:
  using base_t = base_impl_t;
  
public:

  // for standard policy, I need to have and pass the model
  ExplicitStepper(model_t & model,
		  state_t const & y0,
		  res_t const & r0)
    : base_impl_t(model, policy_, y0, r0){}

  // for arbitrary policy, with also the model passed
  ExplicitStepper(model_t & model,
		  pol_t & policyObj,
		  state_t const & y0,
		  res_t const & r0)
    : base_impl_t(model, policyObj, y0, r0){}
  
  ExplicitStepper() = delete;
  ~ExplicitStepper() = default;

private:
  // not used if policy is passed from outside
  pol_std_t policy_; 

};//end class

}} // end namespace rompp::ode
#endif 
