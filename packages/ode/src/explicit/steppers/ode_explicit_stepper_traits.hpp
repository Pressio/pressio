
#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_TRAITS_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_TRAITS_HPP_

#include "../../ode_forward_declarations.hpp"

namespace rompp{ namespace ode{ namespace details{

/*
 * Euler
 */
template<
  typename ode_state_type,
  typename model_type,
  typename ode_residual_type,
  typename ...Args
  >
struct traits<
  ExplicitStepper<ExplicitEnum::Euler, ode_state_type,
		  model_type, ode_residual_type, Args...>
  >{

  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;
  using order_t = unsigned int;
  static constexpr order_t order_value = 1;

  using state_t	   = ode_state_type;
  using residual_t = ode_residual_type;
  using model_t    = model_type;

  // check if scalar is provided in Args
  using ic0 = ::rompp::mpl::variadic::find_if_unary_pred_t<
    std::is_floating_point, Args...>;
  // store the type
  using scalar_t = ::rompp::mpl::variadic::at_or_t<
    void, ic0::value, Args...>;
  static_assert( std::is_void<scalar_t>::value == false,
		 "You need a scalar_type in the ExplicitStepper templates");

  // this is the standard residual policy (just typedef, it is only used
  // if the user does not pass a user-defined policy)
  using standard_res_policy_t = policy::ExplicitResidualStandardPolicy<
    state_t, model_t, residual_t>;

  // check Args if a user-defined residual policy is passed
  using ic1 = ::rompp::mpl::variadic::find_if_unary_pred_t<
    ::rompp::ode::meta::is_legitimate_explicit_residual_policy, Args...>;
  // store the type
  using residual_policy_t = ::rompp::mpl::variadic::at_or_t
    <standard_res_policy_t, ic1::value, Args...>;

  // check if user passed an ops
  using ic2 = ::rompp::mpl::variadic::find_if_quaternary_pred_t<
    scalar_t, state_t, residual_t,
    ::rompp::ode::meta::is_valid_user_defined_ops_for_explicit_ode, Args...>;
  // store the type
  using ops_t = ::rompp::mpl::variadic::at_or_t<void, ic2::value, Args...>;

  using impl_t = impl::ExplicitEulerStepperImpl
    <scalar_t, state_t, model_t, residual_t, residual_policy_t, ops_t>;
};


/*
 * RK4
 */
template<
  typename ode_state_type,
  typename model_type,
  typename ode_residual_type,
  typename ...Args
  >
struct traits<
  ExplicitStepper<ExplicitEnum::RungeKutta4, ode_state_type,
		  model_type, ode_residual_type, Args...>
  >{

  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;
  using order_t = unsigned int;
  static constexpr order_t order_value = 4;

  using state_t	   = ode_state_type;
  using residual_t = ode_residual_type;
  using model_t    = model_type;

  // check if scalar is provided in Args
  using ic0 = ::rompp::mpl::variadic::find_if_unary_pred_t<
    std::is_floating_point, Args...>;
  // store the type
  using scalar_t = ::rompp::mpl::variadic::at_or_t<
    void, ic0::value, Args...>;
  static_assert( std::is_void<scalar_t>::value == false,
		 "You need a scalar_type in the ExplicitStepper templates");

  // this is the standard residual policy (just typedef, it is only used
  // if the user does not pass a user-defined policy)
  using standard_res_policy_t = policy::ExplicitResidualStandardPolicy<
    state_t, model_t, residual_t>;

  // check Args if a user-defined residual policy is passed
  using ic1 = ::rompp::mpl::variadic::find_if_unary_pred_t<
    ::rompp::ode::meta::is_legitimate_explicit_residual_policy, Args...>;
  // store the type
  using residual_policy_t = ::rompp::mpl::variadic::at_or_t
    <standard_res_policy_t, ic1::value, Args...>;

  using impl_t = impl::ExplicitRungeKutta4StepperImpl
    <scalar_t, state_t, model_t, residual_t, residual_policy_t>;
};


}}}//end namespace rompp::ode::details
#endif
