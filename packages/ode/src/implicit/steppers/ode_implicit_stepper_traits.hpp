
#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_TRAITS_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_TRAITS_HPP_

#include "../../ode_forward_declarations.hpp"

namespace rompp{ namespace ode{ namespace details{

template<
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename model_type,
  typename ...Args
  >
struct traits<
  ImplicitStepper<
    ImplicitEnum::Euler,
    state_type, residual_type,
    jacobian_type, model_type,
    Args...>
  > {

  static constexpr ode::ImplicitEnum enum_id = ode::ImplicitEnum::Euler;
  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;

  using state_t		  = state_type;
  using residual_t	  = residual_type;
  using jacobian_t	  = jacobian_type;
  using model_t		  = model_type;
  using aux_stepper_t	  = void;
  using scalar_t	  = typename core::details::traits<state_type>::scalar_t;

  static constexpr unsigned int order_value = 1;
  static constexpr unsigned int steps = 1;

  // standard policies (only used if not passed a user-defined policy)
  using standard_res_policy_t = policy::ImplicitResidualStandardPolicy<
    state_t, model_t, residual_t>;
  using standard_jac_policy_t = policy::ImplicitJacobianStandardPolicy<
    state_t, model_t, jacobian_t>;

  // check Args if a user-defined admissible residual policy is passed
  using ic1 = ::rompp::ode::meta::find_legitimate_implicit_residual_policy_t<
    ImplicitEnum::Euler, 1, state_t, residual_t, model_t, scalar_t, Args...>;
  using residual_policy_t = ::rompp::mpl::variadic::at_or_t
    <standard_res_policy_t, ic1::value, Args...>;

  // check Args if a user-defined admissible jacobian policy is passed
  using ic2 = ::rompp::ode::meta::find_legitimate_implicit_jacobian_policy_t<
    ImplicitEnum::Euler, state_t, jacobian_t, model_t, scalar_t, Args...>;
  using jacobian_policy_t = ::rompp::mpl::variadic::at_or_t
    <standard_jac_policy_t, ic2::value, Args...>;

  using impl_t = impl::ImplicitEulerStepperImpl
    <state_type, residual_type, jacobian_type,
     model_t, residual_policy_t, jacobian_policy_t>;
};


template<
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename model_type,
  typename ...Args
  >
struct traits<
  ImplicitStepper<
    ImplicitEnum::BDF2,
    state_type, residual_type,
    jacobian_type, model_type,
    Args...>
  > {

  static constexpr ode::ImplicitEnum enum_id = ode::ImplicitEnum::BDF2;
  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;

  using this_t = ImplicitStepper<ImplicitEnum::BDF2,
				 state_type, residual_type,
				 jacobian_type, model_type,
				 Args...>;

  using state_t		  = state_type;
  using residual_t	  = residual_type;
  using jacobian_t	  = jacobian_type;
  using model_t		  = model_type;
  using scalar_t	  = typename core::details::traits<state_type>::scalar_t;

  static constexpr unsigned int order_value = 2;
  static constexpr unsigned int steps = 2;

  // for BDF2 the user has to pass an auxiliary stepper
  using ic0 = ::rompp::mpl::variadic::find_if_binary_pred_t<
    this_t, ::rompp::ode::meta::is_legitimate_auxiliary_stepper, Args...>;
  using aux_stepper_t = ::rompp::mpl::variadic::at_or_t<void, ic0::value, Args...>;

  // standard policies (only used if user-defined policies not passed)
  using standard_res_policy_t = policy::ImplicitResidualStandardPolicy<
    state_t, model_t, residual_t>;
  using standard_jac_policy_t = policy::ImplicitJacobianStandardPolicy<
    state_t, model_t, jacobian_t>;

  // check Args if a user-defined admissible residual policy is passed
  using ic1 = ::rompp::ode::meta::find_legitimate_implicit_residual_policy_t<
    ImplicitEnum::BDF2, 2, state_t, residual_t, model_t, scalar_t, Args...>;
  using residual_policy_t = ::rompp::mpl::variadic::at_or_t
    <standard_res_policy_t, ic1::value, Args...>;

  // check Args if a user-defined admissible jacobian policy is passed
  using ic2 = ::rompp::ode::meta::find_legitimate_implicit_jacobian_policy_t<
    ImplicitEnum::BDF2, state_t, jacobian_t, model_t, scalar_t, Args...>;
  using jacobian_policy_t = ::rompp::mpl::variadic::at_or_t
    <standard_jac_policy_t, ic2::value, Args...>;

  using impl_t = impl::ImplicitBDF2StepperImpl
    <state_type, residual_type, jacobian_type,
     model_t, aux_stepper_t, residual_policy_t, jacobian_policy_t>;
};

}}}//end namespace rompp::ode::details
#endif
