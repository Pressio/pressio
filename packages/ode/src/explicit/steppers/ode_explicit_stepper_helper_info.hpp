
#ifndef ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_HELPER_INFO_HPP_
#define ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_HELPER_INFO_HPP_

#include "../policies/ode_explicit_residual_standard_policy.hpp"
#include "../../ode_basic_meta.hpp"
#include "../policies/ode_explicit_policies_meta.hpp"

namespace rompp{ namespace ode{ namespace impl{
  

template <ExplicitSteppersEnum, typename...>
struct stepper_helper_info;
//--------------------------------------------

template <typename... Args>
struct stepper_helper_info<ExplicitSteppersEnum::Euler, Args...>{
  
  // find if a state type is passed
  using ic1_t = ::tinympl::variadic::find_if_t<
    ode::meta::is_legitimate_explicit_state_type, Args...>;
  static constexpr auto k1 = ic1_t::value;
  using state_type = ::tinympl::variadic::at_or_t<void, k1, Args...>;
  static_assert( !std::is_void<state_type>::value,
		 "Invalid state type, it has to be a template");
  
  // find if Args contains a model type
  using ic2_t = ::tinympl::variadic::find_if_t<
    ode::meta::is_legitimate_model_for_explicit_ode, Args...>;
  static constexpr auto k2 = ic2_t::value;
  using model_type = ::tinympl::variadic::at_or_t<void, k2, Args...>;

  // find if a residual type is passed (usually = to state_type)
  using ic3_t = ::tinympl::variadic::find_if_t<
    ode::meta::is_legitimate_explicit_residual_type, Args...>;
  static constexpr auto k3 = ic3_t::value;
  using res_type = ::tinympl::variadic::at_or_t<state_type,k3,Args...>;
  
  // the standard policy is 
  using res_std_pol_type = policy::ExplicitResidualStandardPolicy<
    state_type, model_type, res_type>;
  
  // find if Args contains a residual policy type
  using ic4_t = ::tinympl::variadic::find_if_t<
    ode::meta::is_legitimate_explicit_residual_policy, Args...>;
  static constexpr auto k4 = ic4_t::value;
  using residual_policy_type =
    ::tinympl::variadic::at_or_t<res_std_pol_type, k4, Args...>;
  static_assert( std::is_same<residual_policy_type,
		 res_std_pol_type>::value,"");
  
  // the impl class finally is
  using base_impl_type = impl::ExplicitEulerStepperImpl<
    state_type, model_type, res_type, residual_policy_type>;

};
//--------------------------------------------



// template <ExplicitSteppersEnum,
// 	  typename state_t, typename model_t,
// 	  typename residual_t,
// 	  typename res_policy_t = void,
// 	  typename enable = void>
// struct stepper_helper_info;
// //--------------------------------------------
      
// template <typename state_t,
// 	  typename model_t,
// 	  typename residual_t>
// struct stepper_helper_info<ExplicitSteppersEnum::Euler,
// 		    state_t, model_t, residual_t, void>{
  
//   using res_pol_t =
//     policy::ExplicitResidualStandardPolicy<
//     state_t,model_t,residual_t>;

//   using base_std_t =
//     impl::ExplicitEulerStepperImpl<
//     state_t, model_t, residual_t, res_pol_t>;
// };
// //--------------------------------------------
      
// template <typename state_t,
// 	  typename model_t,
// 	  typename residual_t,
// 	  typename res_policy_t>
// struct stepper_helper_info<ExplicitSteppersEnum::Euler,
// 			   state_t, model_t,
// 			   residual_t, res_policy_t>{
  
//   using res_pol_t = res_policy_t;
  
//   using base_std_t = impl::ExplicitEulerStepperImpl<
//     state_t, model_t, residual_t, res_pol_t>;
// };
// //--------------------------------------------
      
}}} // end namespace rompp::ode::impl
#endif
