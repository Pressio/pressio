
#ifndef ODE_STEPPERS_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_TRAITS_HPP_
#define ODE_STEPPERS_EXPLICIT_STEPPERS_EXPLICIT_STEPPER_TRAITS_HPP_

#include "ode_forward_declarations.hpp"

namespace ode{
namespace details{

template<typename state_type,
	 typename space_residual_type,
	 typename scalar_type,
	 typename model_type,
	 typename residual_policy_type
	 >
struct traits<impl::ExplicitEulerStepperImpl<state_type,
					     space_residual_type,
					     scalar_type,
					     model_type,
					     residual_policy_type>
	      >
{
  using state_t =  state_type;
  using space_residual_t = space_residual_type;
  using scalar_t = scalar_type;
  using model_t = model_type;
  using residual_policy_t = residual_policy_type;

  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;

  using order_t = unsigned int;
  static constexpr order_t order_value = 1;    
};


////////////////////////////////////////////////////////////////


template<typename state_type,
	 typename space_residual_type,
	 typename scalar_type,
	 typename model_type,
	 typename residual_policy_type,
	 typename butcher_table_type>
struct traits<impl::ExplicitRungeKutta4StepperImpl<state_type,
						   space_residual_type,
						   scalar_type,
						   model_type,
						   residual_policy_type,
						   butcher_table_type>
	      >
{
  using state_t =  state_type;
  using space_residual_t = space_residual_type;
  using scalar_t = scalar_type;
  using model_t = model_type;
  using residual_policy_t = residual_policy_type;
  using butcher_table_t = butcher_table_type;
  
  using order_t = unsigned int;
  static constexpr order_t order_value = 4;

  static constexpr bool is_implicit = false;
  static constexpr bool is_explicit = true;
};

  
}//end namespace details
}//end namespace ode
#endif
