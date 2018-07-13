
#ifndef ODE_EXPLICIT_STEPPER_TRAITS_HPP_
#define ODE_EXPLICIT_STEPPER_TRAITS_HPP_

#include "ode_forward_declarations.hpp"

namespace ode{
namespace details{

template<typename state_type,
	 typename residual_type,
	 typename scalar_type,
	 typename model_type,
	 typename time_type,
	 typename sizer_type,
	 typename residual_policy_type
	 >
struct traits<impl::explicitEulerStepperImpl<state_type,
					     residual_type,
					     scalar_type,
					     model_type,
					     time_type,
					     sizer_type,
					     residual_policy_type>
	      >
{
  using state_t =  state_type;
  using residual_t = residual_type;
  using scalar_t = scalar_type;
  using model_t = model_type;
  using time_t = time_type;
  using sizer_t = sizer_type;
  using residual_policy_t = residual_policy_type;
  
  using order_t = unsigned int;
  static constexpr order_t order_value = 1;    
};


////////////////////////////////////////////////////////////////


template<typename state_type,
	 typename residual_type,
	 typename scalar_type,
	 typename model_type,
	 typename time_type,
	 typename sizer_type,
	 typename residual_policy_type,
	 typename butcher_table_type>
struct traits<impl::explicitRungeKutta4StepperImpl<state_type,
						   residual_type,
						   scalar_type,
						   model_type,
						   time_type,
						   sizer_type,
						   residual_policy_type,
						   butcher_table_type>
	      >
{
  using state_t =  state_type;
  using residual_t = residual_type;
  using scalar_t = scalar_type;
  using model_t = model_type;
  using time_t = time_type;
  using sizer_t = sizer_type;
  using residual_policy_t = residual_policy_type;
  using butcher_table_t = butcher_table_type;
  
  using order_t = unsigned int;
  static constexpr order_t order_value = 4;
};

  
}//end namespace details
}//end namespace ode
#endif
