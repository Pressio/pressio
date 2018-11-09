
#ifndef ODE_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_HELPER_INFO_HPP_
#define ODE_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_HELPER_INFO_HPP_

#include "../../ode_basic_meta.hpp"
#include "../../ode_is_legitimate_model_for_implicit_ode.hpp"

#include "../policies/meta/ode_is_legitimate_implicit_jacobian_policy.hpp"
#include "../policies/meta/ode_is_legitimate_implicit_residual_policy.hpp"

#include "../policies/standard/ode_implicit_euler_residual_standard_policy.hpp"
#include "../policies/standard/ode_implicit_euler_jacobian_standard_policy.hpp"
#include "../policies/standard/ode_implicit_bdf2_residual_standard_policy.hpp"
#include "../policies/standard/ode_implicit_bdf2_jacobian_standard_policy.hpp"
#include "../../../../core/src/meta/tinympl/at.hpp"
#include "../../../../core/src/meta/tinympl/find_if.hpp"


namespace rompp{ namespace ode{ namespace impl{

//!!!!!!!!!!!!!!!!!
#ifdef HAVE_CPP14
//!!!!!!!!!!!!!!!!!
      
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

	// find if a jacobian type is passed, this has to be passed
	using ic4_t = ::tinympl::variadic::find_if_t<
	  meta::is_legitimate_jacobian_type, Args...>;
	static constexpr auto k4 = ic4_t::value;
	using jac_type = ::tinympl::variadic::at_or_t<void,k4,Args...>;
	static_assert( !std::is_void<jac_type>::value,
		       "cannot pass void jacob type to implicit stepper, it has to be a valid template");

	// find if a residual type is passed: if not, default = state_type
	// remove the state_type from Args... so that we can find the residual type
	using newargs_t =
	  typename ::tinympl::variadic::erase<k1,k1+1, std::tuple, Args...>::type;  
	using ic3_t = typename ::tinympl::find_if<newargs_t,
						  meta::is_legitimate_implicit_residual_type>::type;
	static constexpr auto k3 = ic3_t::value;
	using res_type = typename ::tinympl::at_or<state_type,k3,newargs_t>::type;

      };
      //--------------------------------------------

      template <ImplicitSteppersEnum, typename...>
      struct implicit_stepper_helper_info;

      //--------------------------------------------

      template <typename... Args>
      struct implicit_stepper_helper_info<ImplicitSteppersEnum::Euler, Args...>
	: implicit_stepper_helper_info_base<Args...>{

	// static constexpr ImplicitSteppersEnum myenum = ImplicitSteppersEnum::Euler;  
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

	// static constexpr ImplicitSteppersEnum myenum = ImplicitSteppersEnum::BDF2;  
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

      
//!!!!!!!!!!!!!!!!!      
#else
//!!!!!!!!!!!!!!!!!

      
      template<ImplicitSteppersEnum whichone,
	       typename ode_state_type,
	       typename ode_residual_type,
	       typename ode_jacobian_type,
	       typename model_type,
	       typename aux_stepper_type,
	       typename residual_policy_type,
	       typename jacobian_policy_type>
      struct implicit_stepper_helper_info;

      //--------------------------------------------

      // backward Euler, standard policies
      template <typename ode_state_type,
		typename ode_residual_type,
		typename ode_jacobian_type,
		typename model_type,
		typename aux_stepper_type>
      struct implicit_stepper_helper_info<ImplicitSteppersEnum::Euler,
					  ode_state_type, ode_residual_type,
					  ode_jacobian_type, model_type,
					  aux_stepper_type, void, void>{
	
	using res_std_pol_type = policy::ImplicitEulerResidualStandardPolicy<
	  ode_state_type, model_type, ode_residual_type>;
	using jac_std_pol_type = policy::ImplicitEulerJacobianStandardPolicy<
	  ode_state_type, model_type, ode_jacobian_type>;
	using base_impl_type = impl::ImplicitEulerStepperImpl<
	  ode_state_type, ode_residual_type, ode_jacobian_type,
	  model_type, res_std_pol_type, jac_std_pol_type>;
      };
      //--------------------------------------------      

      // backward Euler, user-defined policies
      template <typename ode_state_type,
		typename ode_residual_type,
		typename ode_jacobian_type,
		typename model_type,
		typename aux_stepper_type,
		typename residual_policy_type,
		typename jacobian_policy_type>
      struct implicit_stepper_helper_info<ImplicitSteppersEnum::Euler,
					  ode_state_type, ode_residual_type,
					  ode_jacobian_type, model_type,
					  aux_stepper_type, residual_policy_type,
					  jacobian_policy_type>{
	using base_impl_type = impl::ImplicitEulerStepperImpl<
	  ode_state_type, ode_residual_type, ode_jacobian_type,
	  model_type, residual_policy_type, jacobian_policy_type>;
      };
      //--------------------------------------------      

      // BDF2, standard policies
      template <typename ode_state_type,
		typename ode_residual_type,
		typename ode_jacobian_type,
		typename model_type,
		typename aux_stepper_type>
      struct implicit_stepper_helper_info<ImplicitSteppersEnum::BDF2,
					  ode_state_type, ode_residual_type,
					  ode_jacobian_type, model_type,
					  aux_stepper_type, void, void>{
	
	using res_std_pol_type = policy::ImplicitBDF2ResidualStandardPolicy<
	  ode_state_type, model_type, ode_residual_type>;
	using jac_std_pol_type = policy::ImplicitBDF2JacobianStandardPolicy<
	  ode_state_type, model_type, ode_jacobian_type>;
	using base_impl_type = impl::ImplicitBDF2StepperImpl<
	  ode_state_type, ode_residual_type, ode_jacobian_type,
	  model_type, aux_stepper_type, res_std_pol_type, jac_std_pol_type>;
      };
      //--------------------------------------------

      // BDF2, user-defined policies
      template <typename ode_state_type,
		typename ode_residual_type,
		typename ode_jacobian_type,
		typename model_type,
		typename aux_stepper_type,
		typename residual_policy_type,
		typename jacobian_policy_type>
      struct implicit_stepper_helper_info<ImplicitSteppersEnum::BDF2,
					  ode_state_type, ode_residual_type,
					  ode_jacobian_type, model_type,
					  aux_stepper_type, residual_policy_type,
					  jacobian_policy_type>{
	using base_impl_type = impl::ImplicitBDF2StepperImpl<
	  ode_state_type, ode_residual_type, ode_jacobian_type,
	  model_type, aux_stepper_type, residual_policy_type, jacobian_policy_type>;
      };
      //--------------------------------------------      

      
//!!!!!!!!!!!!!!!!!
#endif
//!!!!!!!!!!!!!!!!!      
      
}}} // end namespace rompp::impl
#endif
