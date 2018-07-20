
#ifndef ODE_STEPPERS_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_TRAITS_HPP_
#define ODE_STEPPERS_IMPLICIT_STEPPERS_IMPLICIT_STEPPER_TRAITS_HPP_

#include "ode_forward_declarations.hpp"

namespace ode{
namespace details{
  
template<typename state_type,
	 typename residual_type,
	 typename jacobian_type,
	 typename model_type,
	 typename sizer_type,
	 typename residual_policy_type,
	 typename jacobian_policy_type>
struct traits< impl::ImplicitEulerStepperImpl<
		 state_type, residual_type, jacobian_type,
		 typename core::details::traits<state_type>::scalar_t,
		 model_type, sizer_type, residual_policy_type,
		 jacobian_policy_type>>
{

  using stepper_t =
    impl::ImplicitEulerStepperImpl<state_type,
				   residual_type,
				   jacobian_type,
	     typename core::details::traits<state_type>::scalar_t,
				   model_type,
				   sizer_type,
				   residual_policy_type,
				   jacobian_policy_type>;
  using state_t =  state_type;
  using residual_t = residual_type;
  using jacobian_t =  jacobian_type;
  using scalar_t = typename core::details::traits<state_type>::scalar_t;
  using model_t = model_type;
  using sizer_t = sizer_type;
  using residual_policy_t = residual_policy_type;
  using jacobian_policy_t = jacobian_policy_type;

  static constexpr bool is_implicit = true;
  static constexpr bool is_explicit = false;
  
  using order_t = unsigned int;
  static constexpr order_t order_value = 1;
  static constexpr order_t steps = 1;
};

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

  
}//end namespace details
}//end namespace ode

#endif











// template<typename state_type,
// 	 typename residual_type,
// 	 typename jacobian_type,
// 	 typename scalar_type,
// 	 typename model_type,
// 	 typename time_type,
// 	 typename sizer_type,
// 	 typename solver_policy_type,
// 	 typename aux_start_stepper_type,
// 	 typename residual_policy_type,
// 	 typename jacobian_policy_type>
// struct traits< impl::implicitBDF2StepperImpl<state_type,
// 					     residual_type,
// 					     jacobian_type,
// 					     scalar_type,
// 					     model_type,
// 					     time_type,
// 					     sizer_type,
// 					     solver_policy_type,
// 					     aux_start_stepper_type,
// 					     residual_policy_type,
// 					     jacobian_policy_type> >
// {
//   using stepper_t =
//     impl::implicitBDF2StepperImpl<state_type,
// 				  residual_type,
// 				  jacobian_type,
// 				  scalar_type,
// 				  model_type,
// 				  time_type,
// 				  sizer_type,
// 				  solver_policy_type,
// 				  aux_start_stepper_type,
// 				  residual_policy_type,
// 				  jacobian_policy_type>;
//   using state_t =  state_type;
//   using residual_t = residual_type;
//   using jacobian_t =  jacobian_type;
//   using scalar_t = scalar_type;    
//   using model_t = model_type;
//   using time_t = time_type;
//   using sizer_t = sizer_type;
//   using solver_policy_t = solver_policy_type;
//   using residual_policy_t = residual_policy_type;
//   using jacobian_policy_t = jacobian_policy_type;

//   static constexpr bool advanceIncrement = residual_policy_t::advanceIncrement;
//   static_assert(residual_policy_t::advanceIncrement ==
//   		jacobian_policy_t::advanceIncrement,
//   		"Residual and jacobian policies BOTH need to advance full state or just increment wrt initial condition. In this case they are not");

//   using order_t = unsigned int;
//   static constexpr order_t order_value = 2;
//   static constexpr order_t steps = 2;
// };
  


// ////////////////////////////////////////////////////////
// ////////////////////////////////////////////////////////


// template<typename state_type,
// 	 typename residual_type,
// 	 typename jacobian_type,
// 	 typename scalar_type,
// 	 typename model_type,
// 	 typename time_type,
// 	 typename sizer_type,
// 	 typename solver_policy_type,
// 	 typename aux_start_stepper_type,
// 	 typename residual_policy_type,
// 	 typename jacobian_policy_type>
// struct traits< impl::implicitBDF3StepperImpl<state_type,
// 					     residual_type,
// 					     jacobian_type,
// 					     scalar_type,
// 					     model_type,
// 					     time_type,
// 					     sizer_type,
// 					     solver_policy_type,
// 					     aux_start_stepper_type,
// 					     residual_policy_type,
// 					     jacobian_policy_type> >
// {
//   using stepper_t =
//     impl::implicitBDF3StepperImpl<state_type,
// 				  residual_type,
// 				  jacobian_type,
// 				  scalar_type,
// 				  model_type,
// 				  time_type,
// 				  sizer_type,
// 				  solver_policy_type,
// 				  aux_start_stepper_type,
// 				  residual_policy_type,
// 				  jacobian_policy_type>;
//   using state_t =  state_type;
//   using residual_t = residual_type;
//   using jacobian_t =  jacobian_type;
//   using scalar_t = scalar_type;    
//   using model_t = model_type;
//   using time_t = time_type;
//   using sizer_t = sizer_type;
//   using solver_policy_t = solver_policy_type;
//   using residual_policy_t = residual_policy_type;
//   using jacobian_policy_t = jacobian_policy_type;

//   static constexpr bool advanceIncrement = residual_policy_t::advanceIncrement;
//   static_assert(residual_policy_t::advanceIncrement ==
//   		jacobian_policy_t::advanceIncrement,
//   		"Residual and jacobian policies BOTH need to advance full state or just increment wrt initial condition. In this case they are not");

//   using order_t = unsigned int;
//   static constexpr order_t order_value = 3;
//   static constexpr order_t steps = 3;
// };  
  



// template<typename state_type,
// 	 typename residual_type,
// 	 typename jacobian_type,
// 	 typename scalar_type,
// 	 typename model_type,
// 	 typename time_type,
// 	 typename sizer_type,
// 	 typename solver_policy_type,
// 	 typename residual_policy_type,
// 	 typename jacobian_policy_type>
// struct traits< impl::implicitAdamsMoulton1StepperImpl<state_type,
// 						      residual_type,
// 						      jacobian_type,
// 						      scalar_type,
// 						      model_type,
// 						      time_type,
// 						      sizer_type,
// 						      solver_policy_type,
// 						      residual_policy_type,
// 						      jacobian_policy_type>>
// {
//   using stepper_t =
//     impl::implicitAdamsMoulton1StepperImpl<state_type,
// 					   residual_type,
// 					   jacobian_type,
// 					   scalar_type,
// 					   model_type,
// 					   time_type,
// 					   sizer_type,
// 					   solver_policy_type,
// 					   residual_policy_type,
// 					   jacobian_policy_type>;
//   using state_t =  state_type;
//   using residual_t = residual_type;
//   using jacobian_t =  jacobian_type;
//   using scalar_t = scalar_type;    
//   using model_t = model_type;
//   using time_t = time_type;
//   using sizer_t = sizer_type;
//   using solver_policy_t = solver_policy_type;
//   using residual_policy_t = residual_policy_type;
//   using jacobian_policy_t = jacobian_policy_type;

//   static constexpr bool advanceIncrement = residual_policy_t::advanceIncrement;
//   static_assert(residual_policy_t::advanceIncrement ==
//   		jacobian_policy_t::advanceIncrement,
//   		"Residual and jacobian policies BOTH need to advance full state or just increment wrt initial condition. In this case they are not");

//   using order_t = unsigned int;
//   static constexpr order_t order_value = 2;
//   static constexpr order_t steps = 1;
// };
