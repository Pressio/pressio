
#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_TRAITS_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_TRAITS_HPP_

#include "../../ode_forward_declarations.hpp"

namespace rompp{ namespace ode{ namespace details{


// backward Euler, STANDARD policies
template<
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename model_type
  >
struct traits<
  ImplicitStepper<ImplicitEnum::Euler,
		  state_type, residual_type, jacobian_type,
		  model_type, void, void, void>
  >{

  static constexpr ::rompp::ode::ImplicitEnum enum_id = ImplicitEnum::Euler;
  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;
  using order_t	= unsigned int;
  static constexpr order_t order_value = 1;
  static constexpr order_t steps = 1;

  using state_t		  = state_type;
  using residual_t	  = residual_type;
  using jacobian_t	  = jacobian_type;
  using model_t		  = model_type;
  using aux_stepper_t	  = void;
  using residual_policy_t = policy::ImplicitResidualStandardPolicy<
				state_type, model_type, residual_type>;
  using jacobian_policy_t = policy::ImplicitJacobianStandardPolicy<
				state_type, model_type, jacobian_type>;
  using scalar_t	  = typename core::details::traits<state_type>::scalar_t;

  using impl_t = impl::ImplicitEulerStepperImpl<state_type, residual_type,
						jacobian_type, model_type,
						residual_policy_t,
						jacobian_policy_t>;
};



// backward Euler, user-defined policies
template<
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename model_type,
  typename residual_policy_type,
  typename jacobian_policy_type
  >
struct traits<
  ImplicitStepper<ImplicitEnum::Euler,
		  state_type, residual_type, jacobian_type,
		  model_type, void, residual_policy_type,
		  jacobian_policy_type>
  >{

  static constexpr ::rompp::ode::ImplicitEnum enum_id = ImplicitEnum::Euler;
  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;
  using order_t	= unsigned int;
  static constexpr order_t order_value = 1;
  static constexpr order_t steps = 1;

  using state_t		  = state_type;
  using residual_t	  = residual_type;
  using jacobian_t	  = jacobian_type;
  using model_t		  = model_type;
  using aux_stepper_t	  = void;
  using residual_policy_t = residual_policy_type;
  using jacobian_policy_t = jacobian_policy_type;
  using scalar_t	  = typename core::details::traits<state_type>::scalar_t;

  using impl_t = impl::ImplicitEulerStepperImpl<state_type, residual_type,
						jacobian_type, model_type,
						residual_policy_t,
						jacobian_policy_t>;
};



// BDF2, STANDARD policies
template<
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename model_type,
  typename aux_stepper_type
  >
struct traits<
  ImplicitStepper<ImplicitEnum::BDF2,
		  state_type, residual_type, jacobian_type,
		  model_type, aux_stepper_type, void, void>
  >{

  static constexpr ::rompp::ode::ImplicitEnum enum_id = ImplicitEnum::BDF2;
  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;
  using order_t	= unsigned int;
  static constexpr order_t order_value = 2;
  static constexpr order_t steps = 2;

  using state_t		  = state_type;
  using residual_t	  = residual_type;
  using jacobian_t	  = jacobian_type;
  using model_t		  = model_type;
  using aux_stepper_t	  = aux_stepper_type;
  using residual_policy_t = policy::ImplicitResidualStandardPolicy<
				state_type, model_type, residual_type>;
  using jacobian_policy_t = policy::ImplicitJacobianStandardPolicy<
				state_type, model_type, jacobian_type>;
  using scalar_t	  = typename core::details::traits<state_type>::scalar_t;

  using impl_t = impl::ImplicitBDF2StepperImpl<state_type, residual_type,
  					       jacobian_type, model_type,
  					       aux_stepper_t,
  					       residual_policy_t,
  					       jacobian_policy_t>;
};


// BFD2 Euler, user-defined policies
template<
  typename state_type,
  typename residual_type,
  typename jacobian_type,
  typename model_type,
  typename aux_stepper_type,
  typename residual_policy_type,
  typename jacobian_policy_type
  >
struct traits<
  ImplicitStepper<ImplicitEnum::BDF2,
		  state_type, residual_type, jacobian_type,
		  model_type, aux_stepper_type, residual_policy_type,
		  jacobian_policy_type>
  >{

  static constexpr ::rompp::ode::ImplicitEnum enum_id = ImplicitEnum::BDF2;
  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;
  using order_t	= unsigned int;
  static constexpr order_t order_value = 2;
  static constexpr order_t steps = 2;

  using state_t		  = state_type;
  using residual_t	  = residual_type;
  using jacobian_t	  = jacobian_type;
  using model_t		  = model_type;
  using aux_stepper_t	  = aux_stepper_type;
  using residual_policy_t = residual_policy_type;
  using jacobian_policy_t = jacobian_policy_type;
  using scalar_t	  = typename core::details::traits<state_type>::scalar_t;

  using impl_t = impl::ImplicitBDF2StepperImpl<state_type, residual_type,
					       jacobian_type, model_type,
					       aux_stepper_t,
					       residual_policy_t,
					       jacobian_policy_t>;
};


}}}//end namespace rompp::ode::details
#endif
