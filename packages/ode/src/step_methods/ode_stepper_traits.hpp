
#ifndef ODE_STEPPER_TRAITS_HPP_
#define ODE_STEPPER_TRAITS_HPP_

#include "ode_forward_declarations.hpp"

namespace ode{
namespace details{

template<typename state_type,
	 typename residual_type,
	 typename scalar_type,
	 typename resizer_fctor_type,
	 typename model_type,
	 typename time_type,
	 typename residual_policy_type
	 >
struct traits<explicitEulerStepper<state_type,residual_type,
				   scalar_type,resizer_fctor_type,
				   model_type,time_type,
				   residual_policy_type>
	      >
{
  using stepper_t =
    explicitEulerStepper<state_type,residual_type,
			 scalar_type, resizer_fctor_type,
			 model_type, time_type,
			 residual_policy_type>;
  using state_t =  state_type;
  using residual_t = residual_type;
  using scalar_t = scalar_type;
  using resizer_t = resizer_fctor_type;
  using model_t = model_type;
  using time_t = time_type;
  using residual_policy_t = residual_policy_type;

  using order_t = unsigned int;
  static constexpr order_t order_value = 1;    
};


template<typename state_type,
	 typename residual_type,
	 typename scalar_type,
	 typename resizer_fctor_type,
	 typename model_type,
	 typename time_type,
	 typename residual_policy_type
	 >
struct traits<explicitRungeKutta4Stepper<state_type,residual_type,scalar_type,
					 resizer_fctor_type,model_type,time_type,
					 residual_policy_type> >
{
  using stepper_t = explicitRungeKutta4Stepper<state_type,residual_type,scalar_type,
					       resizer_fctor_type,model_type,time_type,
					       residual_policy_type>;
  using state_t =  state_type;
  using residual_t = residual_type;
  using scalar_t = scalar_type;
  using resizer_t = resizer_fctor_type;
  using model_t = model_type;
  using time_t = time_type;
  using residual_policy_t = residual_policy_type;
  using order_t = unsigned int;
  static constexpr order_t order_value = 4;
};

  
  
template<typename state_type,
	 typename residual_type,
	 typename jacobian_type,
	 typename scalar_type,
	 typename resizer_fctor_type,
	 typename model_type,
	 typename time_type,
	 typename solver_policy_type,
	 typename residual_policy_type,
	 typename jacobian_policy_type>
struct traits< implicitEulerStepper<state_type, residual_type, jacobian_type, scalar_type,
				    resizer_fctor_type, model_type, time_type,
				    solver_policy_type,residual_policy_type,
				    jacobian_policy_type>>
{
  using stepper_t = implicitEulerStepper<state_type, residual_type, jacobian_type,
					 scalar_type, resizer_fctor_type, model_type,
					 time_type,solver_policy_type,residual_policy_type,
					 jacobian_policy_type>;
  using state_t =  state_type;
  using residual_t = residual_type;
  using jacobian_t =  jacobian_type;
  using scalar_t = scalar_type;    
  using resizer_t = resizer_fctor_type;
  using model_t = model_type;
  using time_t = time_type;
  using solver_policy_t = solver_policy_type;
  using residual_policy_t = residual_policy_type;
  using jacobian_policy_t = jacobian_policy_type;

  using order_t = unsigned int;
  static constexpr order_t order_value = 1;
};

  
}//end namespace details
}//end namespace ode

#endif

