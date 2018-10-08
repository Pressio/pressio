
#ifndef ODE_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_HPP_
#define ODE_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_HPP_

#include "ode_implicit_stepper_helper_info.hpp"

namespace rompp{ namespace ode{ 
      
template<ImplicitSteppersEnum whichone, typename... Args>
class ImplicitStepper
  : public impl::implicit_stepper_helper_info<whichone,
					      Args...>::base_impl_type{

  using info_t = impl::implicit_stepper_helper_info<whichone, Args...>;
  using state_t = typename info_t::state_type;
  using model_t = typename info_t::model_type;
  using res_t = typename info_t::res_type;
  using jac_t = typename info_t::jac_type;
  using residual_pol_t = typename info_t::residual_policy_type;
  using residual_pol_std_t = typename info_t::res_std_pol_type;
  using jacobian_pol_t = typename info_t::jacobian_policy_type;
  using jacobian_pol_std_t = typename info_t::jac_std_pol_type;
  using base_impl_t = typename info_t::base_impl_type;

//this needs to be public, it is detected by integrators
public:
  using base_t = base_impl_t;
  
public:
  // for standard policy, I need to have and pass the model
  ImplicitStepper(const model_t & model,
		  state_t const & y0,
		  res_t const & r0)
    : base_impl_t(model,
		  res_policy_obj_,
		  jac_policy_obj_,
		  y0, r0){}

  // when using arbitrary policy, bu also passing the model 
  ImplicitStepper(const model_t & model,
  		  const residual_pol_t & resPolicyObj,
		  const jacobian_pol_t & jacPolicyObj,
  		  state_t const & y0,
  		  res_t const & r0)
    : base_impl_t(model,
		  resPolicyObj,
		  jacPolicyObj,
		  y0, r0){}
  
  ImplicitStepper() = delete;
  ~ImplicitStepper() = default;

private:
  // not used if policy is passed from outside
  residual_pol_t res_policy_obj_;
  jacobian_pol_t jac_policy_obj_;

};//end class

}} // end namespace rompp::ode
#endif 
