
#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_TRAITS_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_TRAITS_HPP_

#include "../../ode_forward_declarations.hpp"

namespace rompp{ namespace ode{ namespace details{

/*
 * Eurler, standard policy
 */
template<
  typename ode_state_type,
  typename model_type,
  typename ode_residual_type
  >
struct traits<
  ExplicitStepper<ExplicitEnum::Euler, ode_state_type,
		  model_type, ode_residual_type, void>
  >{

  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;
  using order_t = unsigned int;
  static constexpr order_t order_value = 1;

  using state_t	   = ode_state_type;
  using residual_t = ode_residual_type;
  using scalar_t   = typename core::details::traits<ode_state_type>::scalar_t;
  using model_t    = model_type;
  using residual_policy_t = policy::ExplicitResidualStandardPolicy<
    ode_state_type, model_type, ode_residual_type>;
  using impl_t = impl::ExplicitEulerStepperImpl<state_t, model_t,
					  residual_t, residual_policy_t>;
};


/*
 * Eurler, user-define policy
 */
template<
  typename ode_state_type,
  typename model_type,
  typename ode_residual_type,
  typename residual_policy_type
  >
struct traits<
  ExplicitStepper<ExplicitEnum::Euler, ode_state_type,
		  model_type, ode_residual_type,
		  residual_policy_type>
  >{

  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;
  using order_t = unsigned int;
  static constexpr order_t order_value = 1;

  using state_t	   = ode_state_type;
  using residual_t = ode_residual_type;
  using scalar_t   = typename core::details::traits<ode_state_type>::scalar_t;
  using model_t    = model_type;
  using residual_policy_t = residual_policy_type;
  using impl_t = impl::ExplicitEulerStepperImpl<state_t, model_t,
					  residual_t, residual_policy_t>;
};



/*
 * RK4, standard policy
 */
template<
  typename ode_state_type,
  typename model_type,
  typename ode_residual_type
  >
struct traits<
  ExplicitStepper<ExplicitEnum::RungeKutta4, ode_state_type,
		  model_type, ode_residual_type, void>
  >{

  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;
  using order_t = unsigned int;
  static constexpr order_t order_value = 4;

  using state_t	   = ode_state_type;
  using residual_t = ode_residual_type;
  using scalar_t   = typename core::details::traits<ode_state_type>::scalar_t;
  using model_t    = model_type;
  using residual_policy_t = policy::ExplicitResidualStandardPolicy<
    ode_state_type, model_type, ode_residual_type>;
  using impl_t = impl::ExplicitRungeKutta4StepperImpl<state_t, model_t,
						residual_t,
						residual_policy_t>;
};


/*
 * RK4, user-defined policy
 */
template<
  typename ode_state_type,
  typename model_type,
  typename ode_residual_type,
  typename residual_policy_type
  >
struct traits<
  ExplicitStepper<ExplicitEnum::RungeKutta4, ode_state_type,
		  model_type, ode_residual_type,
		  residual_policy_type>
  >{

  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;
  using order_t = unsigned int;
  static constexpr order_t order_value = 4;

  using state_t	   = ode_state_type;
  using residual_t = ode_residual_type;
  using scalar_t   = typename core::details::traits<ode_state_type>::scalar_t;
  using model_t    = model_type;
  using residual_policy_t = residual_policy_type;
  using impl_t = impl::ExplicitRungeKutta4StepperImpl<state_t, model_t,
						residual_t,
						residual_policy_t>;
};


}}}//end namespace rompp::ode::details
#endif
