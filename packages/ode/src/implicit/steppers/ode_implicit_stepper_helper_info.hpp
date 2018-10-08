
#ifndef ODE_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_HELPER_INFO_HPP_
#define ODE_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_HELPER_INFO_HPP_

#include "../../ode_basic_meta.hpp"
#include "../policies/meta/ode_implicit_policies_meta.hpp"
#include "../policies/meta/ode_implicit_euler_policies_meta.hpp"
#include "../policies/standard/ode_implicit_euler_residual_standard_policy.hpp"
#include "../policies/standard/ode_implicit_euler_jacobian_standard_policy.hpp"
#include "../policies/meta/ode_implicit_bdf2_policies_meta.hpp"
#include "../policies/standard/ode_implicit_bdf2_residual_standard_policy.hpp"
#include "../policies/standard/ode_implicit_bdf2_jacobian_standard_policy.hpp"


namespace rompp{ namespace ode{ namespace impl{
  
template <typename... Args>
struct implicit_stepper_helper_info_base{
  
  // find if a state type is passed
  using ic1_t = ::tinympl::variadic::find_if_t<
    meta::is_legitimate_implicit_state_type, Args...>;
  static constexpr auto k1 = ic1_t::value;
  using state_type = ::tinympl::variadic::at_or_t<void, k1, Args...>;
  static_assert( !std::is_void<state_type>::value,
  "cannot pass void state type to implicit stepper, it has to be a valid template");
  
  // find if Args contains a model type
  using ic2_t = ::tinympl::variadic::find_if_t<
    meta::is_legitimate_model_for_implicit_ode, Args...>;
  static constexpr auto k2 = ic2_t::value;
  using model_type = ::tinympl::variadic::at_or_t<void, k2, Args...>;

  // find if a residual type is passed: if not, default = state_type
  using ic3_t = ::tinympl::variadic::find_if_t<
    meta::is_legitimate_implicit_residual_type, Args...>;
  static constexpr auto k3 = ic3_t::value;
  using res_type = ::tinympl::variadic::at_or_t<state_type,k3,Args...>;

  // find if a jacobian type is passed, this has to be passed
  using ic4_t = ::tinympl::variadic::find_if_t<
    meta::is_legitimate_jacobian_type, Args...>;
  static constexpr auto k4 = ic4_t::value;
  using jac_type = ::tinympl::variadic::at_or_t<void,k4,Args...>;
  static_assert( !std::is_void<jac_type>::value,
  "cannot pass void jacob type to implicit stepper, it has to be a valid template");
  
};
//--------------------------------------------


template <ImplicitSteppersEnum, typename...>
struct implicit_stepper_helper_info;

//--------------------------------------------

template <typename... Args>
struct implicit_stepper_helper_info<ImplicitSteppersEnum::Euler, Args...>
  : implicit_stepper_helper_info_base<Args...>{

  auto myenum = ImplicitSteppersEnum::Euler;  
  using base_t = implicit_stepper_helper_info_base<Args...>;
  using state_type = typename base_t::state_type;
  using model_type = typename base_t::model_type;
  using res_type = typename base_t::res_type;
  using jac_type = typename base_t::jac_type;
  using auxiliary_stepper_type = void;
  
  // the RESIDUAL standard policy is 
  using res_std_pol_type = policy::ImplicitEulerResidualStandardPolicy<
    state_type, model_type, res_type>;

  // the JACOBIAN standard policy is 
  using jac_std_pol_type = policy::ImplicitEulerJacobianStandardPolicy<
    state_type, model_type, jac_type>;
  
  // find if Args contains a non-standard RESIDUAL policy type
  using ic1_t = ::tinympl::variadic::find_if_t<
    meta::is_legitimate_implicit_euler_residual_policy, Args...>;
  static constexpr auto k1 = ic1_t::value;
  using residual_policy_type =
    ::tinympl::variadic::at_or_t<res_std_pol_type, k1, Args...>;

  // find if Args contains a non-standard JACOBIAN policy type
  using ic2_t = ::tinympl::variadic::find_if_t<
    meta::is_legitimate_implicit_euler_jacobian_policy, Args...>;
  static constexpr auto k2 = ic2_t::value;
  using jacobian_policy_type =
    ::tinympl::variadic::at_or_t<jac_std_pol_type, k2, Args...>;
  
  // the impl class finally is
  using base_impl_type = impl::ImplicitEulerStepperImpl<
    state_type, res_type, jac_type, model_type,
    residual_policy_type, jacobian_policy_type>;

};
//--------------------------------------------      

template <typename... Args>
struct implicit_stepper_helper_info<ImplicitSteppersEnum::BDF2, Args...>
  : implicit_stepper_helper_info_base<Args...>{

  using base_t = implicit_stepper_helper_info_base<Args...>;
  using state_type = typename base_t::state_type;
  using model_type = typename base_t::model_type;
  using res_type = typename base_t::res_type;
  using jac_type = typename base_t::jac_type;
  
  // the RESIDUAL standard policy is 
  using res_std_pol_type = policy::ImplicitBDF2ResidualStandardPolicy<
    state_type, model_type, res_type>;

  // the JACOBIAN standard policy is 
  using jac_std_pol_type = policy::ImplicitBDF2JacobianStandardPolicy<
    state_type, model_type, jac_type>;
  
  // find if Args contains a non-standard RESIDUAL policy type
  using ic1_t = ::tinympl::variadic::find_if_t<
    meta::is_legitimate_implicit_bdf2_residual_policy, Args...>;
  static constexpr auto k1 = ic1_t::value;
  using residual_policy_type =
    ::tinympl::variadic::at_or_t<res_std_pol_type, k1, Args...>;

  // find if Args contains a non-standard JACOBIAN policy type
  using ic2_t = ::tinympl::variadic::find_if_t<
    meta::is_legitimate_implicit_bdf2_jacobian_policy, Args...>;
  static constexpr auto k2 = ic2_t::value;
  using jacobian_policy_type =
    ::tinympl::variadic::at_or_t<jac_std_pol_type, k2, Args...>;

  // find if Args contains an auxiliary stepper for first few steps
  // maybe set one by default
  using ic3_t = ::tinympl::variadic::find_if_t<
    meta::is_legitimate_auxiliary_stepper, Args...>;
  static constexpr auto k3 = ic3_t::value;
  using auxiliary_stepper_type =
    ::tinympl::variadic::at_or_t<void, k3, Args...>;
  static_assert( !std::is_void<auxiliary_stepper_type>::value,
  "cannot pass void aux stepper type to BDF2, it has to be a valid template");
  
  // the impl class finally is
  using base_impl_type = impl::ImplicitBDF2StepperImpl<
    state_type, res_type, jac_type, model_type, 
    auxiliary_stepper_type, residual_policy_type, 
    jacobian_policy_type>;

};
//--------------------------------------------      

      
      
}}} // end namespace rompp::impl
#endif
