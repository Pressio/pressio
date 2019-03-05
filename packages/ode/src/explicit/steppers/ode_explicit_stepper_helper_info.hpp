
#ifndef ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_HELPER_INFO_HPP_
#define ODE_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_HELPER_INFO_HPP_

#include "../../meta/ode_basic_meta.hpp"
#include "../../meta/ode_is_legitimate_model_for_explicit_ode.hpp"
#include "../policies/ode_explicit_policies_meta.hpp"
#include "../policies/ode_explicit_residual_standard_policy.hpp"

namespace rompp{ namespace ode{ namespace impl{

//!!!!!!!!!!!!!!!!!
#ifdef HAVE_CPP14
//!!!!!!!!!!!!!!!!!
      /*
       * TODO: once we move to c++14, all this needs
       * to be transfereed/adapted to the traits class
       */

      template <typename... Args>
      struct explicit_stepper_helper_info_base{

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

      };
      //--------------------------------------------


      template <ExplicitEnum, typename...>
      struct explicit_stepper_helper_info;
      //--------------------------------------------

      template <typename... Args>
      struct explicit_stepper_helper_info<ExplicitEnum::Euler, Args...>
	: explicit_stepper_helper_info_base<Args...>{

	using base_t = explicit_stepper_helper_info_base<Args...>;
	using state_type = typename base_t::state_type;
	using model_type = typename base_t::model_type;
	using res_type = typename base_t::res_type;
	using res_std_pol_type = typename base_t::res_std_pol_type;
	using residual_policy_type = typename base_t::residual_policy_type;

	// the impl class finally is
	using base_impl_type = impl::ExplicitEulerStepperImpl<
	  state_type, model_type, res_type, residual_policy_type>;

      };
      //--------------------------------------------

      template <typename... Args>
      struct explicit_stepper_helper_info<ExplicitEnum::RungeKutta4, Args...>
	: explicit_stepper_helper_info_base<Args...>{

	using base_t = explicit_stepper_helper_info_base<Args...>;
	using state_type = typename base_t::state_type;
	using model_type = typename base_t::model_type;
	using res_type = typename base_t::res_type;
	using res_std_pol_type = typename base_t::res_std_pol_type;
	using residual_policy_type = typename base_t::residual_policy_type;

	// the impl class finally is
	using base_impl_type = impl::ExplicitRungeKutta4StepperImpl<
	  state_type, model_type, res_type, residual_policy_type>;

      };

//!!!!!!!!!!!!!!!!
#endif //HAVE_CPP14
//!!!!!!!!!!!!!!!!

}}} // end namespace rompp::ode::impl
#endif
